use crate::{
    aggregation::compute_aggregation,
    enrich::enrichment_testing,
    io::{GeneFrame, SgrnaFrame, SimpleFrame},
    model::model_mean_variance,
    norm::normalize_counts,
    utils::{config::Configuration, logging::Logger},
};
use anyhow::Result;

/// Performs the `MAGeCK` Differential Expression and Gene Aggregation Algorithm
pub fn mageck(
    frame: &SimpleFrame,
    labels_controls: &[String],
    labels_treatments: &[String],
    config: &Configuration,
    logger: &Logger,
) -> Result<()> {
    let labels = [labels_controls, labels_treatments].concat();
    let count_matrix = frame.data_matrix(&labels)?;
    let sgrna_names = frame.get_sgrna_names();
    let gene_names = frame.get_gene_names();
    let n_controls = labels_controls.len();

    logger.start_mageck();
    logger.num_sgrnas(sgrna_names);
    logger.num_genes(gene_names);
    logger.norm_method(config.normalization());
    logger.aggregation_method(config.aggregation());
    logger.correction(config.correction());

    // Normalize
    let normed_matrix = normalize_counts(&count_matrix, config.normalization(), logger);

    // Mean-Variance Modeling
    let adj_var = model_mean_variance(&normed_matrix, n_controls, config.model_choice(), logger);

    // sgRNA Ranking (Enrichment)
    let sgrna_results =
        enrichment_testing(&normed_matrix, &adj_var, n_controls, config.correction());

    // Gene Ranking (Aggregation)
    let aggregation_results = compute_aggregation(
        config.aggregation(),
        &sgrna_results,
        gene_names,
        logger,
        config.correction(),
    );

    // Build sgRNA DataFrame
    let sgrna_frame = SgrnaFrame::new(sgrna_names, gene_names, &adj_var, &sgrna_results);
    sgrna_frame.write(config.prefix())?;

    // Build Gene DataFrame
    let gene_frame = GeneFrame::new(&aggregation_results);
    gene_frame.write(config.prefix())?;

    Ok(())
}
