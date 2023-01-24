use std::fs::File;
use polars::prelude::{CsvReader, CsvWriter, SerReader, SerWriter, DataFrame, PolarsError, df, Series, NamedFrom};
use ndarray::Array1;
use crate::{
    enrich::EnrichmentResult,
    aggregation::AggregationResult
};

/// Reads a dataframe from file
pub fn load_dataframe(path: &str) -> Result<DataFrame, PolarsError>
{
    CsvReader::from_path(path)?
        .has_header(true)
        .with_delimiter(b'\t')
        .finish()
}

/// Writes a dataframe to file
fn write_dataframe(name: &str, frame: &mut DataFrame) -> Result<(), PolarsError>
{
    let file = File::create(name)?;
    CsvWriter::new(file)
        .has_header(true)
        .with_delimiter(b'\t')
        .finish(frame)
}

/// Write the `sgRNA` results dataframe given a prefix
pub fn write_sgrna_results(prefix: &str, frame: &mut DataFrame) -> Result<(), PolarsError>
{
    let filename = format!("{}.sgrna_results.tab", prefix);
    write_dataframe(&filename, frame)
}

/// Write the Gene results dataframe given a prefix
pub fn write_gene_results(prefix: &str, frame: &mut DataFrame) -> Result<(), PolarsError>
{
    let filename = format!("{}.gene_results.tab", prefix);
    write_dataframe(&filename, frame)
}

/// Builds the gene-level aggregation results dataframe
pub fn build_gene_dataframe(aggregation_results: &AggregationResult) -> Result<DataFrame, PolarsError>
{
    df!(
        "gene" => aggregation_results.genes(),
        "score_low" => aggregation_results.score_low().to_vec(),
        "pvalue_low" => aggregation_results.pvalues_low().to_vec(),
        "fdr_low" => aggregation_results.fdr_low().to_vec(),
        "score_high" => aggregation_results.score_high().to_vec(),
        "pvalue_high" => aggregation_results.pvalues_high().to_vec(),
        "fdr_high" => aggregation_results.fdr_high().to_vec(),
    )
}

/// Builds the sgrna-level differential expression results dataframe
pub fn build_sgrna_dataframe(
    sgrna_names: &[String],
    gene_names: &[String],
    adj_var: &Array1<f64>,
    sgrna_results: &EnrichmentResult) -> Result<DataFrame, PolarsError>
{
    df!(
        "sgrna" => &sgrna_names,
        "gene" => gene_names,
        "control" => sgrna_results.control_means().to_vec(),
        "treatment" => sgrna_results.treatment_means().to_vec(),
        "adj_var" => adj_var.to_vec(),
        "pvalue_low" => sgrna_results.pvalues_low().to_vec(),
        "pvalue_high" => sgrna_results.pvalues_high().to_vec(),
        "pvalue_twosided" => sgrna_results.pvalues_twosided().to_vec(),
        "fdr" => sgrna_results.fdr().to_vec()
    )
}

