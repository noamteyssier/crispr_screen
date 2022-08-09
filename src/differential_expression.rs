use anyhow::Result;
use std::fs::File;
use ndarray::{s, Array1};
use polars::prelude::{DataFrame, Series, NamedFrom, df, CsvWriter, SerWriter};
use crate::{
    utils::{
        parse_to_string_vec, parse_to_ndarray, 
        model_mean_variance, enrichment_testing},
    rra::alpha_rra,
    norm::{Normalization, normalize_counts}
};


pub enum GeneAggregation {
    AlpaRRA{alpha: f64, npermutations: usize}
}

fn compute_aggregation(
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

fn write_table(
    name: &str,
    frame: &mut DataFrame) -> Result<()>
{
    let file = File::create(name)?;
    CsvWriter::new(file)
        .has_header(true)
        .with_delimiter(b'\t')
        .finish(frame)?;
    Ok(())
}

fn write_sgrna_results(
    prefix: &str,
    frame: &mut DataFrame) -> Result<()>
{
    let filename = format!("{}.sgrna_results.tab", prefix);
    write_table(&filename, frame)
}

fn write_gene_results(
    prefix: &str,
    frame: &mut DataFrame) -> Result<()>
{
    let filename = format!("{}.gene_results.tab", prefix);
    write_table(&filename, frame)
}

pub fn mageck(
    frame: &DataFrame,
    labels_controls: &[String],
    labels_treatments: &[String],
    prefix: &str,
    normalization: Normalization,
    aggregation: GeneAggregation,
    ) -> Result<()>
{
    let columns = frame.get_column_names();
    let labels = [labels_controls, labels_treatments].concat();
    let count_matrix = parse_to_ndarray(frame, &labels)?;
    let sgrna_names = parse_to_string_vec(frame, columns[0])?;
    let gene_names = parse_to_string_vec(frame, columns[1])?;

    // Normalize
    let normed_matrix = normalize_counts(&count_matrix, normalization);

    // Mean-Variance Modeling
    let adj_var = model_mean_variance(&normed_matrix, labels_controls.len());

    // sgRNA Ranking (Enrichment)
    let (sgrna_pvalues_low, sgrna_pvalues_high)= enrichment_testing(&normed_matrix, &adj_var, labels_controls.len());

    // Gene Ranking (Aggregation)
    let (genes, gene_pvalues_low, gene_pvalues_high) = compute_aggregation(
        aggregation,
        &sgrna_pvalues_low,
        &sgrna_pvalues_high,
        &gene_names);

    let mut sgrna_frame = df!(
        "sgrna" => &sgrna_names,
        "gene" => gene_names,
        "control" => normed_matrix.slice(s![.., 0]).to_vec(),
        "treatment" => normed_matrix.slice(s![.., 1]).to_vec(),
        "adj_var" => adj_var.to_vec(),
        "pvalues_low" => sgrna_pvalues_low.to_vec(),
        "pvalues_high" => sgrna_pvalues_high.to_vec()
    )?;

    let mut gene_frame = df!(
        "gene" => genes,
        "pvalues_low" => gene_pvalues_low.to_vec(),
        "pvalues_high" => gene_pvalues_high.to_vec()
    )?;

    write_sgrna_results(prefix, &mut sgrna_frame)?;
    write_gene_results(prefix, &mut gene_frame)?;


    Ok(())
}
