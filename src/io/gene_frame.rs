use std::{fs::File, io::{BufWriter, Write}};
use anyhow::Result;
use ndarray::Array1;
use crate::aggregation::AggregationResult;


pub struct GeneFrame<'a> {
    gene: &'a [String],
    score_low: &'a Array1<f64>,
    pvalue_low: &'a Array1<f64>,
    fdr_low: &'a Array1<f64>,
    score_high: &'a Array1<f64>,
    pvalue_high: &'a Array1<f64>,
    fdr_high: &'a Array1<f64>,
    size: usize,
}

impl <'a> GeneFrame <'a> {

    pub fn new(aggregation_results: &'a AggregationResult) -> Self {
        Self {
            gene: aggregation_results.genes(),
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
        let mut writer = File::create(format!("{}.gene_results.tab", prefix))
            .map(BufWriter::new)?;

        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            "gene",
            "score_low",
            "pvalue_low",
            "fdr_low",
            "score_high",
            "pvalue_high",
            "fdr_high",
        )?;

        for idx in 0..self.size {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.gene[idx],
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
