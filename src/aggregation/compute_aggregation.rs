use ndarray::Array1;
use super::alpha_rra;

pub enum GeneAggregation {
    AlpaRRA{alpha: f64, npermutations: usize}
}

pub fn compute_aggregation(
    agg: GeneAggregation,
    sgrna_pvalues_low: &Array1<f64>,
    sgrna_pvalues_high: &Array1<f64>,
    gene_names: &Vec<String>) -> (Vec<String>, Array1<f64>, Array1<f64>)
{
    match agg {
        GeneAggregation::AlpaRRA { alpha, npermutations } => {
            let (genes, gene_pvalues_low) = alpha_rra(&sgrna_pvalues_low, &gene_names, alpha, npermutations);
            let (_, gene_pvalues_high) = alpha_rra(&sgrna_pvalues_high, &gene_names, alpha, npermutations);
            (genes, gene_pvalues_low, gene_pvalues_high)
        }
    }
}
