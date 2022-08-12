use ndarray::Array1;
use crate::utils::logging::Logger;

use super::{alpha_rra, inc};

/// Enum describing the different gene aggregation procedures and their associated configurations.
#[derive(Debug)]
pub enum GeneAggregation <'a> {
    AlpaRRA{ alpha: f64, npermutations: usize },
    Inc { token: &'a str }
}

/// Computes gene aggregation using the provided method and associated configurations.
pub fn compute_aggregation(
    agg: &GeneAggregation,
    sgrna_pvalues_low: &Array1<f64>,
    sgrna_pvalues_high: &Array1<f64>,
    gene_names: &Vec<String>,
    logger: &Logger) -> (Vec<String>, Array1<f64>, Array1<f64>)
{
    logger.start_gene_aggregation();
    match agg {
        GeneAggregation::AlpaRRA { alpha, npermutations } => {
            let (genes, gene_pvalues_low) = alpha_rra(sgrna_pvalues_low, gene_names, *alpha, *npermutations, logger);
            let (_, gene_pvalues_high) = alpha_rra(sgrna_pvalues_high, gene_names, *alpha, *npermutations, logger);
            (genes, gene_pvalues_low, gene_pvalues_high)
        },
        GeneAggregation::Inc { token } => {
            let (genes, gene_pvalues_low) = inc(sgrna_pvalues_low, gene_names, token, logger);
            let (_, gene_pvalues_high) = inc(sgrna_pvalues_high, gene_names, token, logger);
            (genes, gene_pvalues_low, gene_pvalues_high)
        }
    }
}
