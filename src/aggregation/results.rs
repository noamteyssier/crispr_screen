use adjustp::{adjust, Procedure};
use ndarray::Array1;

pub struct AggregationResult {
    genes: Vec<String>,
    gene_fc: Array1<f64>,
    gene_log2_fc: Array1<f64>,
    pvalues_low: Array1<f64>,
    pvalues_high: Array1<f64>,
    aggregation_score_low: Array1<f64>,
    aggregation_score_high: Array1<f64>,
    fdr_low: Array1<f64>,
    fdr_high: Array1<f64>,
}
impl AggregationResult {
    pub fn new(
        genes: Vec<String>,
        gene_fc: Array1<f64>,
        pvalues_low: Array1<f64>,
        pvalues_high: Array1<f64>,
        aggregation_score_low: Array1<f64>,
        aggregation_score_high: Array1<f64>,
        correction: Procedure,
    ) -> Self {
        let fdr_low = Self::fdr_adjustment(&pvalues_low, correction);
        let fdr_high = Self::fdr_adjustment(&pvalues_high, correction);
        let gene_log2_fc = Self::calculate_log_fold_change(&gene_fc);
        Self {
            genes,
            gene_fc,
            gene_log2_fc,
            pvalues_low,
            pvalues_high,
            aggregation_score_low,
            aggregation_score_high,
            fdr_low,
            fdr_high,
        }
    }

    fn fdr_adjustment(pvalues: &Array1<f64>, correction: Procedure) -> Array1<f64> {
        Array1::from_vec(adjust(pvalues.as_slice().unwrap(), correction))
    }

    fn calculate_log_fold_change(gene_fc: &Array1<f64>) -> Array1<f64> {
        gene_fc.mapv(|x| x.log2())
    }

    pub fn genes(&self) -> &Vec<String> {
        &self.genes
    }

    pub fn gene_fc(&self) -> &Array1<f64> {
        &self.gene_fc
    }

    pub fn gene_log2_fc(&self) -> &Array1<f64> {
        &self.gene_log2_fc
    }

    pub fn pvalues_low(&self) -> &Array1<f64> {
        &self.pvalues_low
    }

    pub fn pvalues_high(&self) -> &Array1<f64> {
        &self.pvalues_high
    }

    pub fn score_low(&self) -> &Array1<f64> {
        &self.aggregation_score_low
    }

    pub fn score_high(&self) -> &Array1<f64> {
        &self.aggregation_score_high
    }

    pub fn fdr_low(&self) -> &Array1<f64> {
        &self.fdr_low
    }

    pub fn fdr_high(&self) -> &Array1<f64> {
        &self.fdr_high
    }
}
