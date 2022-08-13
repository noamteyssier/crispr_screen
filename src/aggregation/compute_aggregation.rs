use ndarray::{Array1, Array2, Axis};
use crate::utils::logging::Logger;
use hashbrown::HashSet;
use super::{alpha_rra, inc, AggregationResult};

/// Enum describing the different gene aggregation procedures and their associated configurations.
#[derive(Debug)]
pub enum GeneAggregation <'a> {
    AlpaRRA{ alpha: f64, npermutations: usize },
    Inc { token: &'a str }
}

/// Return all indices where values are above zero
fn mask_zeros(array: &Array1<f64>, logger: &Logger) -> HashSet<usize>
{
    let mask = array
        .iter()
        .enumerate()
        .filter(|(_idx, x)| **x > 0.)
        .map(|(idx, _)| idx)
        .collect::<HashSet<usize>>();
    logger.num_zeros(array.len() - mask.len());
    mask
}


/// Filter sgRNAs with zero counts in both samples
fn filter_zeros(
    normed_matrix: &Array2<f64>,
    gene_names: &Vec<String>,
    sgrna_pvalues_low: &Array1<f64>,
    sgrna_pvalues_high: &Array1<f64>,
    logger: &Logger) -> (Vec<String>, Array1<f64>, Array1<f64>)
{
    let sgrna_means = normed_matrix.mean_axis(Axis(1)).unwrap();
    let passing_indices = mask_zeros(&sgrna_means, logger);
    let mut sorted_indices = Vec::from_iter(passing_indices.iter().cloned());
    sorted_indices.sort_unstable();

    let passing_gene_names = sorted_indices.iter().map(|x| gene_names[*x].clone()).collect::<Vec<String>>();
    let passing_sgrna_pvalues_low = sorted_indices.iter().map(|x| sgrna_pvalues_low[*x]).collect::<Array1<f64>>();
    let passing_sgrna_pvalues_high = sorted_indices.iter().map(|x| sgrna_pvalues_high[*x]).collect::<Array1<f64>>();

    (
        passing_gene_names,
        passing_sgrna_pvalues_low,
        passing_sgrna_pvalues_high
    )
}

/// Computes gene aggregation using the provided method and associated configurations.
pub fn compute_aggregation(
    agg: &GeneAggregation,
    normed_matrix: &Array2<f64>,
    sgrna_pvalues_low: &Array1<f64>,
    sgrna_pvalues_high: &Array1<f64>,
    gene_names: &Vec<String>,
    logger: &Logger) -> AggregationResult
{
    logger.start_gene_aggregation();
    
    let (passing_gene_names, passing_sgrna_pvalues_low, passing_sgrna_pvalues_high) = filter_zeros(
        normed_matrix, gene_names, sgrna_pvalues_low, sgrna_pvalues_high, logger);

    match agg {
        GeneAggregation::AlpaRRA { alpha, npermutations } => {
            let (genes, gene_scores_low, gene_pvalues_low) = alpha_rra(
                &passing_sgrna_pvalues_low, 
                &passing_gene_names, 
                *alpha, 
                *npermutations, 
                logger);
            let (_, gene_scores_high, gene_pvalues_high) = alpha_rra(
                &passing_sgrna_pvalues_high, 
                &passing_gene_names, 
                *alpha, 
                *npermutations, 
                logger);
            AggregationResult::new(genes, gene_pvalues_low, gene_pvalues_high, gene_scores_low, gene_scores_high)
        },
        GeneAggregation::Inc { token } => {
            let (genes, gene_scores_low, gene_pvalues_low) = inc(
                &passing_sgrna_pvalues_low, 
                &passing_gene_names, 
                token, 
                logger);
            let (_, gene_scores_high, gene_pvalues_high) = inc(
                &passing_sgrna_pvalues_high, 
                &passing_gene_names, 
                token, 
                logger);
            AggregationResult::new(genes, gene_pvalues_low, gene_pvalues_high, gene_scores_low, gene_scores_high)
        }
    }
}
