use crate::aggregation::AggregationResult;
use anyhow::Result;
use ndarray::Array1;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

pub struct GeneFrame<'a> {
    gene: &'a [String],
    gene_fc: &'a Array1<f64>,
    gene_log2_fc: &'a Array1<f64>,
    score_low: &'a Array1<f64>,
    pvalue_low: &'a Array1<f64>,
    fdr_low: &'a Array1<f64>,
    score_high: &'a Array1<f64>,
    pvalue_high: &'a Array1<f64>,
    fdr_high: &'a Array1<f64>,
    pvalue: &'a Array1<f64>,
    fdr: &'a Array1<f64>,
    phenotype_score: &'a Array1<f64>,
    size: usize,
}

impl<'a> GeneFrame<'a> {
    pub fn new(aggregation_results: &'a AggregationResult) -> Self {
        Self {
            gene: aggregation_results.genes(),
            gene_fc: aggregation_results.gene_fc(),
            gene_log2_fc: aggregation_results.gene_log2_fc(),
            score_low: aggregation_results.score_low(),
            pvalue_low: aggregation_results.pvalues_low(),
            fdr_low: aggregation_results.fdr_low(),
            score_high: aggregation_results.score_high(),
            pvalue_high: aggregation_results.pvalues_high(),
            fdr_high: aggregation_results.fdr_high(),
            pvalue: aggregation_results.pvalue(),
            fdr: aggregation_results.fdr(),
            size: aggregation_results.genes().len(),
            phenotype_score: aggregation_results.phenotype_score(),
        }
    }

    pub fn write(&self, prefix: &str) -> Result<()> {
        let mut writer = File::create(format!("{prefix}.gene_results.tab")).map(BufWriter::new)?;

        writeln!(
            writer,
            "gene\tfold_change\tlog_fold_change\tscore_low\tpvalue_low\tfdr_low\tscore_high\tpvalue_high\tfdr_high\tpvalue\tfdr\tphenotype_score"
        )?;

        for idx in 0..self.size {
            writeln!(
                writer,
                "{}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}",
                self.gene[idx],
                self.gene_fc[idx],
                self.gene_log2_fc[idx],
                self.score_low[idx],
                self.pvalue_low[idx],
                self.fdr_low[idx],
                self.score_high[idx],
                self.pvalue_high[idx],
                self.fdr_high[idx],
                self.pvalue[idx],
                self.fdr[idx],
                self.phenotype_score[idx],
            )?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod testing {
    use crate::aggregation::AggregationResult;
    use ndarray::array;

    #[test]
    fn test_gene_frame() {
        let gene = vec!["gene1".to_string(), "gene2".to_string()];
        let gene_fc = array![1.0, 2.0];
        let score_low = array![0.0, 1.0];
        let pvalue_low = array![0.0, 1.0];
        let fdr_low = array![0.0, 1.0];
        let score_high = array![0.0, 1.0];
        let pvalue_high = array![0.0, 1.0];
        let fdr_high = array![0.0, 1.0];
        let threshold_low = 0.0;
        let threshold_high = 1.0;

        let aggregation_results = AggregationResult::new(
            gene,
            gene_fc,
            pvalue_low,
            pvalue_high,
            fdr_low,
            fdr_high,
            score_low,
            score_high,
            Some(threshold_low),
            Some(threshold_high),
        );
        let gene_frame = super::GeneFrame::new(&aggregation_results);
        assert_eq!(gene_frame.size, 2);
        gene_frame.write("test").unwrap();
    }
}
