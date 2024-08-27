use crate::{
    aggregation::compute_aggregation,
    enrich::enrichment_testing,
    io::{
        get_string_column, match_headers_from_regex_set, validate_ntc, write_gene_frame,
        write_hit_list, write_sgrna_dataframe, Screenviz,
    },
    model::model_mean_variance,
    norm::normalize_counts,
    utils::{config::Configuration, filter::filter_low_counts, logging::Logger},
};
use anyhow::Result;
use polars::prelude::*;
use regex::Regex;

/// Performs the `MAGeCK` Differential Expression and Gene Aggregation Algorithm
pub fn mageck(
    frame: &DataFrame,
    regex_controls: &[Regex],
    regex_treatments: &[Regex],
    config: &Configuration,
    logger: &Logger,
    skip_agg: bool,
) -> Result<()> {
    let control_labels = match_headers_from_regex_set(frame, regex_controls)?;
    let treatment_labels = match_headers_from_regex_set(frame, regex_treatments)?;
    let n_controls = control_labels.len();
    let labels = [control_labels.clone(), treatment_labels.clone()].concat();

    let count_matrix = frame
        .select(&labels)?
        .to_ndarray::<Float64Type>(IndexOrder::Fortran)?;
    let sgrna_names = get_string_column(frame, 0);
    let gene_names = get_string_column(frame, 1);
    validate_ntc(&sgrna_names, config.aggregation())?;

    logger.start_mageck();
    logger.group_names(&control_labels, &treatment_labels);
    logger.num_sgrnas(&sgrna_names);
    logger.num_genes(&gene_names);
    logger.norm_method(config.normalization());
    logger.aggregation_method(config.aggregation());
    logger.correction(*config.correction());

    let normed_matrix = normalize_counts(&count_matrix, config.normalization(), logger);

    // Filter Low Counts
    let (filt_matrix, filt_sgrna_names, filt_gene_names) = filter_low_counts(
        &normed_matrix,
        &sgrna_names,
        &gene_names,
        *config.min_base_mean(),
        logger,
    );

    // Mean-Variance Modeling
    let adj_var = model_mean_variance(&filt_matrix, n_controls, config.model_choice(), logger);

    // sgRNA Ranking (Enrichment)
    let sgrna_results = enrichment_testing(
        &filt_matrix,
        &adj_var,
        n_controls,
        *config.correction(),
        *config.strategy(),
        &logger,
    );

    // Write sgRNA DataFrame
    write_sgrna_dataframe(
        &filt_sgrna_names,
        &filt_gene_names,
        adj_var.as_slice().unwrap(),
        &sgrna_results,
        config.prefix(),
    )?;

    if skip_agg {
        Ok(())
    } else {
        // Gene Ranking (Aggregation)
        let aggregation_results = compute_aggregation(
            config.aggregation(),
            &sgrna_results,
            &filt_gene_names,
            logger,
            *config.correction(),
            *config.seed(),
        );

        // Build Gene DataFrame
        write_gene_frame(&aggregation_results, config.prefix())?;

        // Write hit list
        write_hit_list(&aggregation_results, config, logger)?;

        // Write screenviz config
        let screenviz = Screenviz::new(&aggregation_results, config);
        screenviz.write(config.prefix())?;

        Ok(())
    }
}
