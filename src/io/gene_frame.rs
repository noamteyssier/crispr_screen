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
            size: aggregation_results.genes().len(),
        }
    }

    pub fn write(&self, prefix: &str) -> Result<()> {
        let mut writer = File::create(format!("{prefix}.gene_results.tab")).map(BufWriter::new)?;

        writeln!(
            writer,
            "gene\tfold_change\tlog_fold_change\tscore_low\tpvalue_low\tfdr_low\tscore_high\tpvalue_high\tfdr_high",
        )?;

        for idx in 0..self.size {
            writeln!(
                writer,
                "{}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}",
                self.gene[idx],
                self.gene_fc[idx],
                self.gene_log2_fc[idx],
                self.score_low[idx],
                self.pvalue_low[idx],
                self.fdr_low[idx],
                self.score_high[idx],
                self.pvalue_high[idx],
                self.fdr_high[idx],
            )?;
        }

        Ok(())
    }
}
