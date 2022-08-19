use anyhow::Result;
use polars::prelude::DataFrame;
use crate::{
    utils::{
        io::{write_gene_results, write_sgrna_results, build_gene_dataframe, build_sgrna_dataframe}, 
        parse_to_string_vec, parse_to_ndarray, logging::Logger},
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
    logger: &Logger
    ) -> Result<()>
{
    let columns = frame.get_column_names();
    let labels = [labels_controls, labels_treatments].concat();
    let count_matrix = parse_to_ndarray(frame, &labels)?;
    let sgrna_names = parse_to_string_vec(frame, columns[0])?;
    let gene_names = parse_to_string_vec(frame, columns[1])?;

    logger.start_mageck();
    logger.num_sgrnas(&sgrna_names);
    logger.num_genes(&gene_names);
    logger.norm_method(normalization);
    logger.aggregation_method(aggregation);

    // Normalize
    let normed_matrix = normalize_counts(&count_matrix, normalization);

    // Mean-Variance Modeling
    let adj_var = model_mean_variance(&normed_matrix, labels_controls.len(), logger);

    // sgRNA Ranking (Enrichment)
    let sgrna_results = enrichment_testing(&normed_matrix, &adj_var, labels_controls.len());

    // Gene Ranking (Aggregation)
    let aggregation_results = compute_aggregation(
        aggregation,
        &normed_matrix,
        &sgrna_results,
        &gene_names,
        logger);

    // Build sgRNA DataFrame
    let mut sgrna_frame = build_sgrna_dataframe(
        &sgrna_names, 
        &gene_names, 
        &normed_matrix, 
        &adj_var, 
        &sgrna_results)?;
    write_sgrna_results(prefix, &mut sgrna_frame)?;


    // Build Gene DataFrame
    let mut gene_frame = build_gene_dataframe(&aggregation_results)?;
    write_gene_results(prefix, &mut gene_frame)?;

    Ok(())
}
