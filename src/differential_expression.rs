use anyhow::Result;
use ndarray::s;
use polars::prelude::{DataFrame, Series, NamedFrom, df};
use crate::{
    utils::{
        io::{write_gene_results, write_sgrna_results}, 
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
    let (sgrna_pvalues_low, sgrna_pvalues_high)= enrichment_testing(&normed_matrix, &adj_var, labels_controls.len());

    // Gene Ranking (Aggregation)
    let aggregation_results = compute_aggregation(
        aggregation,
        &normed_matrix,
        &sgrna_pvalues_low,
        &sgrna_pvalues_high,
        &gene_names,
        logger);

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
        "gene" => aggregation_results.genes(),
        "score_low" => aggregation_results.score_low().to_vec(),
        "pvalues_low" => aggregation_results.pvalues_low().to_vec(),
        "fdr_low" => aggregation_results.fdr_low().to_vec(),
        "score_high" => aggregation_results.score_high().to_vec(),
        "pvalues_high" => aggregation_results.pvalues_high().to_vec(),
        "fdr_high" => aggregation_results.fdr_high().to_vec(),
    )?;

    write_sgrna_results(prefix, &mut sgrna_frame)?;
    write_gene_results(prefix, &mut gene_frame)?;

    Ok(())
}
