use anyhow::Result;
use std::fs::File;
use ndarray::s;
use polars::prelude::{DataFrame, Series, NamedFrom, df, CsvWriter, SerWriter};
use crate::{
    utils::{
        parse_to_string_vec, parse_to_ndarray, 
        model_mean_variance, enrichment_testing,
        normalize_counts, Normalization},
    rra::alpha_rra
};


pub fn mageck(
    frame: &DataFrame,
    labels_controls: &[String],
    labels_treatments: &[String],
    normalization: Normalization
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

    let (genes_low, gene_pvalues_low) = alpha_rra(&sgrna_pvalues_low, &gene_names, 0.05, 100);
    let (genes_high, gene_pvalues_high) = alpha_rra(&sgrna_pvalues_high, &gene_names, 0.05, 100);

    let mut sgrna_frame = df!(
        "sgrna" => &sgrna_names,
        "gene" => gene_names,
        "control" => normed_matrix.slice(s![.., 0]).to_vec(),
        "treatment" => normed_matrix.slice(s![.., 1]).to_vec(),
        "adj_var" => adj_var.to_vec(),
        "pvalues_low" => sgrna_pvalues_low.to_vec(),
        "pvalues_high" => sgrna_pvalues_high.to_vec()
    )?;

    let file = File::create("sgrna_results.tab")?;
    CsvWriter::new(file)
        .has_header(true)
        .with_delimiter(b'\t')
        .finish(&mut sgrna_frame)?;

    let mut gene_frame = df!(
        "gene" => genes_low,
        "gene_test" => genes_high,
        "pvalues_low" => gene_pvalues_low.to_vec(),
        "pvalues_high" => gene_pvalues_high.to_vec()
    )?;

    let file = File::create("gene_results.tab")?;
    CsvWriter::new(file)
        .has_header(true)
        .with_delimiter(b'\t')
        .finish(&mut gene_frame)?;
    
    // Gene Ranking (Aggregation)
    Ok(())
}
