use adjustp::{adjust, Procedure};
use ndarray::{Array1, Zip};

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
    phenotype_score: Array1<f64>,
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
        let phenotype_score = Self::calculate_phenotype_score(&fdr_low, &fdr_high, &gene_log2_fc);
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
            phenotype_score,
        }
    }

    fn fdr_adjustment(pvalues: &Array1<f64>, correction: Procedure) -> Array1<f64> {
        Array1::from_vec(adjust(pvalues.as_slice().unwrap(), correction))
    }

    fn calculate_log_fold_change(gene_fc: &Array1<f64>) -> Array1<f64> {
        gene_fc.mapv(f64::log2)
    }

    fn calculate_phenotype_score(
        fdr_low: &Array1<f64>,
        fdr_high: &Array1<f64>,
        gene_log2_fc: &Array1<f64>,
    ) -> Array1<f64> {
        Zip::from(fdr_low)
            .and(fdr_high)
            .and(gene_log2_fc)
            .map_collect(|fdr_low, fdr_high, gene_log2_fc| {
                let fdr = fdr_low.min(*fdr_high);
                let nlfdr = -(fdr.log10());
                nlfdr * gene_log2_fc
            })
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

    pub fn phenotype_score(&self) -> &Array1<f64> {
        &self.phenotype_score
    }
}

#[cfg(test)]
mod testing {
    use super::AggregationResult;
    use adjustp::Procedure;
    use ndarray::Array1;

    #[test]
    fn test_aggregation_result() {
        let genes = vec!["gene1".to_string(), "gene2".to_string()];
        let gene_fc = Array1::from(vec![1.0, 2.0]);
        let pvalues_low = Array1::from(vec![0.1, 0.2]);
        let pvalues_high = Array1::from(vec![0.3, 0.4]);
        let aggregation_score_low = Array1::from(vec![0.5, 0.6]);
        let aggregation_score_high = Array1::from(vec![0.7, 0.8]);
        let correction = Procedure::BenjaminiHochberg;
        let result = AggregationResult::new(
            genes,
            gene_fc,
            pvalues_low,
            pvalues_high,
            aggregation_score_low,
            aggregation_score_high,
            correction,
        );
        assert_eq!(
            result.genes(),
            &vec!["gene1".to_string(), "gene2".to_string()]
        );
        assert_eq!(result.gene_fc(), &Array1::from(vec![1.0, 2.0]));
        assert_eq!(result.gene_log2_fc(), &Array1::from(vec![0.0, 1.0]));
        assert_eq!(result.pvalues_low(), &Array1::from(vec![0.1, 0.2]));
        assert_eq!(result.pvalues_high(), &Array1::from(vec![0.3, 0.4]));
        assert_eq!(result.score_low(), &Array1::from(vec![0.5, 0.6]));
        assert_eq!(result.score_high(), &Array1::from(vec![0.7, 0.8]));
        assert_eq!(result.fdr_low(), &Array1::from(vec![0.2, 0.2]));
        assert_eq!(result.fdr_high(), &Array1::from(vec![0.4, 0.4]));
        assert_eq!(
            result.phenotype_score(),
            &Array1::from(vec![0.0, 0.6989700043360187])
        );
    }
}
