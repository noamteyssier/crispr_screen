use super::{alpha_rra, inc, AggregationResult};
use crate::{
    enrich::EnrichmentResult,
    utils::{agg::aggregate_fold_changes, logging::Logger},
};
use adjustp::Procedure;
use hashbrown::HashSet;
use ndarray::{Array1, Array2, Axis};

/// Enum describing the different gene aggregation procedures and their associated configurations.
#[derive(Debug)]
pub enum GeneAggregation<'a> {
    AlpaRRA {
        alpha: f64,
        npermutations: usize,
        adjust_alpha: bool,
    },
    Inc {
        token: &'a str,
    },
}

/// Return all indices where values are above zero
fn mask_zeros(array: &Array1<f64>, logger: &Logger) -> HashSet<usize> {
    let mask = array
        .iter()
        .enumerate()
        .filter(|(_idx, x)| **x > 0.)
        .map(|(idx, _)| idx)
        .collect::<HashSet<usize>>();
    logger.num_zeros(array.len() - mask.len());
    mask
}

/// Filter `sgRNAs` with zero counts in both samples
fn filter_zeros(
    normed_matrix: &Array2<f64>,
    gene_names: &[String],
    sgrna_pvalues_low: &Array1<f64>,
    sgrna_pvalues_high: &Array1<f64>,
    logger: &Logger,
) -> (Vec<String>, Array1<f64>, Array1<f64>) {
    let sgrna_means = normed_matrix.mean_axis(Axis(1)).unwrap();
    let passing_indices = mask_zeros(&sgrna_means, logger);
    let mut sorted_indices = passing_indices.iter().copied().collect::<Vec<usize>>();
    sorted_indices.sort_unstable();

    let passing_gene_names = sorted_indices
        .iter()
        .map(|x| gene_names[*x].clone())
        .collect::<Vec<String>>();
    let passing_sgrna_pvalues_low = sorted_indices
        .iter()
        .map(|x| sgrna_pvalues_low[*x])
        .collect::<Array1<f64>>();
    let passing_sgrna_pvalues_high = sorted_indices
        .iter()
        .map(|x| sgrna_pvalues_high[*x])
        .collect::<Array1<f64>>();

    (
        passing_gene_names,
        passing_sgrna_pvalues_low,
        passing_sgrna_pvalues_high,
    )
}

/// Calculates an empirical alpha threshold for RRA
fn calculate_empirical_alpha(pvalue_arr: &Array1<f64>, alpha: f64) -> f64 {
    pvalue_arr
        .map(|x| if *x < alpha { 1.0 } else { 0.0 })
        .mean()
        .expect("Error calculating mean in empirical alpha")
}

/// Sets the alpha threshold empirically or returns the untouched alpha
fn set_alpha_threshold(
    pvalue_low: &Array1<f64>,
    pvalue_high: &Array1<f64>,
    alpha: f64,
    adjust_alpha: bool,
) -> (f64, f64) {
    if adjust_alpha {
        (
            calculate_empirical_alpha(pvalue_low, alpha),
            calculate_empirical_alpha(pvalue_high, alpha),
        )
    } else {
        (alpha, alpha)
    }
}

/// Computes gene aggregation using the provided method and associated configurations.
pub fn compute_aggregation(
    agg: &GeneAggregation,
    normed_matrix: &Array2<f64>,
    sgrna_results: &EnrichmentResult,
    gene_names: &[String],
    logger: &Logger,
    correction: Procedure,
) -> AggregationResult {
    logger.start_gene_aggregation();

    let gene_fc_hashmap = aggregate_fold_changes(
        gene_names,
        sgrna_results.fold_change(),
        sgrna_results.pvalues_twosided(),
    );

    let (passing_gene_names, passing_sgrna_pvalues_low, passing_sgrna_pvalues_high) = filter_zeros(
        normed_matrix,
        gene_names,
        sgrna_results.pvalues_low(),
        sgrna_results.pvalues_high(),
        logger,
    );

    let (genes, gene_scores_low, gene_pvalues_low, gene_scores_high, gene_pvalues_high) = match agg
    {
        GeneAggregation::AlpaRRA {
            alpha,
            npermutations,
            adjust_alpha,
        } => {
            let (alpha_low, alpha_high) = set_alpha_threshold(
                &passing_sgrna_pvalues_low,
                &passing_sgrna_pvalues_high,
                *alpha,
                *adjust_alpha,
            );
            logger.report_rra_alpha(alpha_low, alpha_high);

            let (genes, gene_scores_low, gene_pvalues_low) = alpha_rra(
                &passing_sgrna_pvalues_low,
                &passing_gene_names,
                alpha_low,
                *npermutations,
                logger,
            );

            let (_, gene_scores_high, gene_pvalues_high) = alpha_rra(
                &passing_sgrna_pvalues_high,
                &passing_gene_names,
                alpha_high,
                *npermutations,
                logger,
            );

            (
                genes,
                gene_scores_low,
                gene_pvalues_low,
                gene_scores_high,
                gene_pvalues_high,
            )
        }
        GeneAggregation::Inc { token } => {
            let (genes, gene_scores_low, gene_pvalues_low) = inc(
                &passing_sgrna_pvalues_low,
                &passing_gene_names,
                token,
                logger,
            );
            let (_, gene_scores_high, gene_pvalues_high) = inc(
                &passing_sgrna_pvalues_high,
                &passing_gene_names,
                token,
                logger,
            );
            (
                genes,
                gene_scores_low,
                gene_pvalues_low,
                gene_scores_high,
                gene_pvalues_high,
            )
        }
    };

    let gene_fc = genes
        .iter()
        .map(|x| gene_fc_hashmap.get(x).unwrap_or(&0.0))
        .copied()
        .collect::<Array1<f64>>();

    AggregationResult::new(
        genes,
        gene_fc,
        gene_pvalues_low,
        gene_pvalues_high,
        gene_scores_low,
        gene_scores_high,
        correction,
    )
}

#[cfg(test)]
mod testing {
    use super::{calculate_empirical_alpha, filter_zeros, mask_zeros};
    use crate::utils::logging::Logger;
    use ndarray::{Array1, Array2, Axis};
    use ndarray_rand::{
        rand_distr::{Binomial, Uniform},
        RandomExt,
    };

    #[test]
    fn test_empirical_alpha() {
        let array = Array1::from_vec(vec![0.1, 0.2, 0.25, 0.5, 0.5]);
        let alpha = 0.3;
        let emp = calculate_empirical_alpha(&array, alpha);
        assert_eq!(emp, 3.0 / 5.0);
    }

    #[test]
    fn test_mask_zeros() {
        let array = Array1::from_vec(vec![0., 1., 0., 1., 0.]);
        let logger = Logger::new();
        let mask = mask_zeros(&array, &logger);
        assert_eq!(mask.len(), 2);
        assert!(mask.contains(&1));
        assert!(mask.contains(&3));
    }

    #[test]
    fn test_filter_zeros() {
        let logger = Logger::new();
        let array = Array2::random((100, 2), Binomial::new(1, 0.2).unwrap()).mapv(|x| x as f64);
        let means = array.mean_axis(Axis(1)).unwrap();
        let nonzero = mask_zeros(&means, &logger);
        let gene_names = (0..100)
            .map(|x| format!("gene_{x}"))
            .collect::<Vec<String>>();
        let p_low = Array1::random(100, Uniform::new(0.0, 1.0));
        let p_high = Array1::random(100, Uniform::new(0.0, 1.0));
        let (pgn, ppl, pph) = filter_zeros(&array, &gene_names, &p_low, &p_high, &logger);

        assert_eq!(pgn.len(), nonzero.len());
        assert_eq!(ppl.len(), nonzero.len());
        assert_eq!(pph.len(), nonzero.len());
    }
}
