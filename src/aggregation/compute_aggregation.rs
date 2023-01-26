use super::{alpha_rra, inc, AggregationResult};
use crate::{enrich::EnrichmentResult, utils::logging::Logger};
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

    AggregationResult::new(
        genes,
        gene_pvalues_low,
        gene_pvalues_high,
        gene_scores_low,
        gene_scores_high,
        correction,
    )
}
