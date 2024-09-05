use bon::bon;
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
    pvalue: Array1<f64>,
    fdr: Array1<f64>,
    phenotype_score: Array1<f64>,
    threshold_low: Option<f64>,
    threshold_high: Option<f64>,
}
#[bon]
impl AggregationResult {
    #[builder]
    pub fn new(
        genes: Vec<String>,
        gene_fc: Array1<f64>,
        pvalues_low: Array1<f64>,
        pvalues_high: Array1<f64>,
        fdr_low: Array1<f64>,
        fdr_high: Array1<f64>,
        aggregation_score_low: Array1<f64>,
        aggregation_score_high: Array1<f64>,
        threshold_low: Option<f64>,
        threshold_high: Option<f64>,
    ) -> Self {
        let pvalue = Self::select_pvalue(&pvalues_low, &pvalues_high);
        let fdr = Self::select_fdr(&fdr_low, &fdr_high);
        let gene_log2_fc = Self::calculate_log_fold_change(&gene_fc);
        let phenotype_score = Self::calculate_phenotype_score(&pvalue, &gene_log2_fc);
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
            pvalue,
            fdr,
            phenotype_score,
            threshold_low,
            threshold_high,
        }
    }

    fn select_fdr(fdr_low: &Array1<f64>, fdr_high: &Array1<f64>) -> Array1<f64> {
        Zip::from(fdr_low)
            .and(fdr_high)
            .map_collect(|fdr_low, fdr_high| fdr_low.min(*fdr_high))
    }

    fn select_pvalue(pvalue_low: &Array1<f64>, pvalue_high: &Array1<f64>) -> Array1<f64> {
        Zip::from(pvalue_low)
            .and(pvalue_high)
            .map_collect(|pvalue_low, pvalue_high| pvalue_low.min(*pvalue_high))
    }

    fn calculate_log_fold_change(gene_fc: &Array1<f64>) -> Array1<f64> {
        gene_fc.mapv(f64::log2)
    }

    fn calculate_phenotype_score(fdr: &Array1<f64>, gene_log2_fc: &Array1<f64>) -> Array1<f64> {
        Zip::from(fdr)
            .and(gene_log2_fc)
            .map_collect(|fdr, gene_log2_fc| {
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

    pub fn pvalue(&self) -> &Array1<f64> {
        &self.pvalue
    }

    pub fn fdr(&self) -> &Array1<f64> {
        &self.fdr
    }

    pub fn phenotype_score(&self) -> &Array1<f64> {
        &self.phenotype_score
    }

    pub fn threshold_low(&self) -> Option<f64> {
        self.threshold_low
    }

    pub fn threshold_high(&self) -> Option<f64> {
        self.threshold_high
    }
}

#[cfg(test)]
mod testing {
    use super::AggregationResult;
    use ndarray::Array1;

    #[test]
    fn test_aggregation_result() {
        let genes = vec!["gene1".to_string(), "gene2".to_string()];
        let gene_fc = Array1::from(vec![1.0, 2.0]);
        let pvalues_low = Array1::from(vec![0.1, 0.2]);
        let pvalues_high = Array1::from(vec![0.3, 0.4]);
        let fdr_low = Array1::from(vec![0.5, 0.6]);
        let fdr_high = Array1::from(vec![0.7, 0.2]);
        let aggregation_score_low = Array1::from(vec![0.5, 0.6]);
        let aggregation_score_high = Array1::from(vec![0.7, 0.8]);
        let result = AggregationResult::builder()
            .genes(genes)
            .gene_fc(gene_fc)
            .pvalues_low(pvalues_low)
            .pvalues_high(pvalues_high)
            .fdr_low(fdr_low)
            .fdr_high(fdr_high)
            .aggregation_score_low(aggregation_score_low)
            .aggregation_score_high(aggregation_score_high)
            .build();

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
        assert_eq!(result.fdr_low(), &Array1::from(vec![0.5, 0.6]));
        assert_eq!(result.fdr_high(), &Array1::from(vec![0.7, 0.2]));
        assert_eq!(result.pvalue(), &Array1::from(vec![0.1, 0.2]));
        assert_eq!(result.fdr(), &Array1::from(vec![0.5, 0.2]));
        assert_eq!(
            result.phenotype_score(),
            &Array1::from(vec![0.0, 0.6989700043360187])
        );
    }
}
