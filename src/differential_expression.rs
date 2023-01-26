use adjustp::Procedure;
use anyhow::Result;
use crate::{
    io::{SimpleFrame, SgrnaFrame, GeneFrame},
    utils::logging::Logger,
    norm::{Normalization, normalize_counts},
    model::{model_mean_variance, ModelChoice},
    enrich::enrichment_testing,
    aggregation::{GeneAggregation, compute_aggregation}
};

/// Performs the `MAGeCK` Differential Expression and Gene Aggregation Algorithm
pub fn mageck(
    frame: &SimpleFrame,
    labels_controls: &[String],
    labels_treatments: &[String],
    prefix: &str,
    normalization: &Normalization,
    aggregation: &GeneAggregation,
    logger: &Logger,
    correction: &Procedure,
    model_choice: &ModelChoice) -> Result<()>
{
    let labels = [labels_controls, labels_treatments].concat();
    let count_matrix = frame.data_matrix(&labels)?;
    let sgrna_names = frame.get_sgrna_names();
    let gene_names = frame.get_gene_names();

    logger.start_mageck();
    logger.num_sgrnas(sgrna_names);
    logger.num_genes(gene_names);
    logger.norm_method(normalization);
    logger.aggregation_method(aggregation);
    logger.correction(correction);

    // Normalize
    let normed_matrix = normalize_counts(&count_matrix, normalization, logger);

    // Mean-Variance Modeling
    let adj_var = model_mean_variance(&normed_matrix, labels_controls.len(), model_choice, logger);

    // sgRNA Ranking (Enrichment)
    let sgrna_results = enrichment_testing(
        &normed_matrix, 
        &adj_var, 
        labels_controls.len(),
        correction);

    // Gene Ranking (Aggregation)
    let aggregation_results = compute_aggregation(
        aggregation,
        &normed_matrix,
        &sgrna_results,
        &gene_names,
        logger,
        correction);

    // Build sgRNA DataFrame
    let sgrna_frame = SgrnaFrame::new(
        sgrna_names, 
        gene_names, 
        &adj_var, 
        &sgrna_results);
    sgrna_frame.write(prefix)?;


    // Build Gene DataFrame
    let gene_frame = GeneFrame::new(&aggregation_results);
    gene_frame.write(prefix)?;

    Ok(())
}
