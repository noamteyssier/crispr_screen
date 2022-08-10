use anyhow::Result;
use ndarray::s;
use polars::prelude::{DataFrame, Series, NamedFrom, df};
use crate::{
    utils::{
        io::{write_gene_results, write_sgrna_results}, 
        parse_to_string_vec, parse_to_ndarray},
    norm::{Normalization, normalize_counts},
    model::model_mean_variance,
    enrich::enrichment_testing,
    aggregation::{GeneAggregation, compute_aggregation}
};

/// Performs the `MAGeCK` Differential Expression and Gene Aggregation Algorithm
pub fn mageck(
    frame: &DataFrame,
    labels_controls: &[String],
    labels_treatments: &[String],
    prefix: &str,
    normalization: &Normalization,
    aggregation: &GeneAggregation,
    ) -> Result<()>
{
    let columns = frame.get_column_names();
    let labels = [labels_controls, labels_treatments].concat();
    let count_matrix = parse_to_ndarray(frame, &labels)?;
    let sgrna_names = parse_to_string_vec(frame, columns[0])?;
    let gene_names = parse_to_string_vec(frame, columns[1])?;

    // Normalize
    let normed_matrix = normalize_counts(&count_matrix, &normalization);

    // Mean-Variance Modeling
    let adj_var = model_mean_variance(&normed_matrix, labels_controls.len());

    // sgRNA Ranking (Enrichment)
    let (sgrna_pvalues_low, sgrna_pvalues_high)= enrichment_testing(&normed_matrix, &adj_var, labels_controls.len());

    // Gene Ranking (Aggregation)
    let (genes, gene_pvalues_low, gene_pvalues_high) = compute_aggregation(
        &aggregation,
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
