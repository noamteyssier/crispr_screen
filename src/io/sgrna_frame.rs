use crate::enrich::EnrichmentResult;
use anyhow::Result;
use ndarray::Array1;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

pub struct SgrnaFrame<'a> {
    sgrna_names: &'a [String],
    gene_names: &'a [String],
    adj_var: &'a Array1<f64>,
    control: &'a Array1<f64>,
    treatment: &'a Array1<f64>,
    pvalue_low: &'a Array1<f64>,
    pvalue_high: &'a Array1<f64>,
    pvalue_twosided: &'a Array1<f64>,
    fdr: &'a Array1<f64>,
    size: usize,
}
impl<'a> SgrnaFrame<'a> {
    pub fn new(
        sgrna_names: &'a [String],
        gene_names: &'a [String],
        adj_var: &'a Array1<f64>,
        sgrna_results: &'a EnrichmentResult,
    ) -> Self {
        Self {
            sgrna_names,
            gene_names,
            adj_var,
            control: sgrna_results.control_means(),
            treatment: sgrna_results.treatment_means(),
            pvalue_low: sgrna_results.pvalues_low(),
            pvalue_high: sgrna_results.pvalues_high(),
            pvalue_twosided: sgrna_results.pvalues_twosided(),
            fdr: sgrna_results.fdr(),
            size: sgrna_names.len(),
        }
    }

    pub fn write(&self, prefix: &str) -> Result<()> {
        let mut writer = File::create(format!("{prefix}.sgrna_results.tab")).map(BufWriter::new)?;

        writeln!(
            writer,
            "sgrna\tgene\tcontrol\ttreatment\tadj_var\tpvalue_low\tpvalue_high\tpvalue_twosided\tfdr",
        )?;

        for idx in 0..self.size {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                self.sgrna_names[idx],
                self.gene_names[idx],
                self.control[idx],
                self.treatment[idx],
                self.adj_var[idx],
                self.pvalue_low[idx],
                self.pvalue_high[idx],
                self.pvalue_twosided[idx],
                self.fdr[idx]
            )?;
        }

        Ok(())
    }
}
