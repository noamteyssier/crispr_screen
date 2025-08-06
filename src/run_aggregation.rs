use anyhow::Result;
use ndarray::s;
use polars::frame::DataFrame;

use crate::{
    aggregation::compute_aggregation,
    cli::SgrnaColumns,
    enrich::EnrichmentResult,
    io::{get_string_column, to_ndarray, write_gene_frame, write_hit_list, Screenviz},
    utils::{config::Configuration, logging::Logger},
};

pub fn run_aggregation(
    frame: &DataFrame,
    columns: SgrnaColumns,
    config: &Configuration,
    logger: &Logger,
) -> Result<()> {
    let sgrna_names = get_string_column(frame, 0);
    let gene_names = get_string_column(frame, 1);
    let sgrna_matrix = to_ndarray(
        frame,
        &[
            columns.pvalue_low,
            columns.pvalue_high,
            columns.control_mean,
            columns.treatment_mean,
        ],
    )?;
    let enrichment_result = EnrichmentResult::new(
        sgrna_matrix.slice(s![.., 0]).to_owned(),
        sgrna_matrix.slice(s![.., 1]).to_owned(),
        sgrna_matrix.slice(s![.., 2]).to_owned(),
        sgrna_matrix.slice(s![.., 3]).to_owned(),
        *config.correction(),
    );

    logger.start_mageck();
    logger.num_sgrnas(&sgrna_names);
    logger.num_genes(&gene_names);
    logger.aggregation_method(config.aggregation());
    logger.correction(*config.correction());

    let aggregation_results = compute_aggregation(
        config.aggregation(),
        &enrichment_result,
        &gene_names,
        logger,
        *config.correction(),
        *config.seed(),
    )?;

    // Write outputs
    write_gene_frame(&aggregation_results, config.prefix())?;

    // Write hit list
    write_hit_list(&aggregation_results, config, logger)?;

    // Write screenviz config
    let screenviz = Screenviz::new(&aggregation_results, config);
    screenviz.write(config.prefix())?;

    Ok(())
}
