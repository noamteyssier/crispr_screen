use anyhow::Result;

use crate::{
    aggregation::compute_aggregation,
    cli::SgrnaColumns,
    enrich::EnrichmentResult,
    io::{GeneFrame, HitList, Screenviz, SimpleFrame},
    utils::{config::Configuration, logging::Logger},
};

pub fn run_aggregation(
    frame: &SimpleFrame,
    columns: SgrnaColumns,
    config: &Configuration,
    logger: &Logger,
) -> Result<()> {
    let pvalues_low = frame.get_f64_column(&columns.pvalue_low)?;
    let pvalues_high = frame.get_f64_column(&columns.pvalue_high)?;
    let control_means = frame.get_f64_column(&columns.control_mean)?;
    let treatment_means = frame.get_f64_column(&columns.treatment_mean)?;
    let gene_names = frame.get_gene_names();
    let enrichment_result = EnrichmentResult::new(
        pvalues_low,
        pvalues_high,
        control_means,
        treatment_means,
        config.correction(),
    );

    let aggregation_results = compute_aggregation(
        config.aggregation(),
        &enrichment_result,
        gene_names,
        logger,
        config.correction(),
        config.seed(),
    );

    // Build Gene DataFrame
    let gene_frame = GeneFrame::new(&aggregation_results);
    gene_frame.write(config.prefix())?;

    // Write hit list
    let hit_list = HitList::new(&aggregation_results, config);
    logger.hit_list(&hit_list);
    hit_list.write(config.prefix())?;

    // Write screenviz config
    let screenviz = Screenviz::new(&aggregation_results, config);
    screenviz.write(config.prefix())?;

    Ok(())
}
