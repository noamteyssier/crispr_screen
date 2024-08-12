use anyhow::Result;
use polars::prelude::*;
use std::{fs::File, io::BufWriter};

use crate::enrich::EnrichmentResult;

fn build_sgrna_dataframe(
    sgrna_names: &[String],
    gene_names: &[String],
    adj_var: &[f64],
    sgrna_results: &EnrichmentResult,
) -> Result<DataFrame, PolarsError> {
    df!(
        "sgrna" => sgrna_names,
        "gene" => gene_names,
        "base" => sgrna_results.base_means().to_vec(),
        "control" => sgrna_results.control_means().to_vec(),
        "treatment" => sgrna_results.treatment_means().to_vec(),
        "adj_var" => adj_var.to_vec(),
        "log2fc" => sgrna_results.log_fold_change().to_vec(),
        "pvalue_low" => sgrna_results.pvalues_low().to_vec(),
        "pvalue_high" => sgrna_results.pvalues_high().to_vec(),
        "pvalue_twosided" => sgrna_results.pvalues_twosided().to_vec(),
        "fdr" => sgrna_results.fdr().to_vec(),
        "product" => sgrna_results.product().to_vec(),
    )
}

pub fn write_sgrna_dataframe(
    sgrna_names: &[String],
    gene_names: &[String],
    adj_var: &[f64],
    sgrna_results: &EnrichmentResult,
    prefix: &str,
) -> Result<(), PolarsError> {
    let mut df = build_sgrna_dataframe(sgrna_names, gene_names, adj_var, sgrna_results)?;
    df.sort_in_place(["fdr"], Default::default())?;
    let writer = File::create(format!("{}.sgrna_results.tsv", prefix)).map(BufWriter::new)?;
    CsvWriter::new(writer)
        .with_separator(b'\t')
        .include_header(true)
        .with_quote_style(QuoteStyle::Never)
        .with_float_scientific(Some(true))
        .finish(&mut df)
}
