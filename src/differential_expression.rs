use crate::{
    aggregation::compute_aggregation,
    enrich::enrichment_testing,
    io::{GeneFrame, HitList, Screenviz, SgrnaFrame, SimpleFrame},
    model::model_mean_variance,
    norm::normalize_counts,
    utils::{config::Configuration, filter::filter_low_counts, logging::Logger},
};
use anyhow::Result;
use regex::Regex;

/// Performs the `MAGeCK` Differential Expression and Gene Aggregation Algorithm
pub fn mageck(
    frame: &SimpleFrame,
    regex_controls: &[Regex],
    regex_treatments: &[Regex],
    config: &Configuration,
    logger: &Logger,
    skip_agg: bool,
) -> Result<()> {
    let control_labels = frame.match_headers_from_regex_set(regex_controls)?;
    let treatment_labels = frame.match_headers_from_regex_set(regex_treatments)?;
    let n_controls = control_labels.len();
    let labels = [control_labels.clone(), treatment_labels.clone()].concat();
    let count_matrix = frame.data_matrix(&labels)?;
    let sgrna_names = frame.get_sgrna_names();
    let gene_names = frame.get_gene_names();

    frame.validate_ntc(config.aggregation())?;
    logger.start_mageck();
    logger.group_names(&control_labels, &treatment_labels);
    logger.num_sgrnas(sgrna_names);
    logger.num_genes(gene_names);
    logger.norm_method(config.normalization());
    logger.aggregation_method(config.aggregation());
    logger.correction(config.correction());

    // Normalize
    let normed_matrix = normalize_counts(&count_matrix, config.normalization(), logger);

    // Filter Low Counts
    let (filt_matrix, filt_sgrna_names, filt_gene_names) = filter_low_counts(
        &normed_matrix,
        sgrna_names,
        gene_names,
        config.min_base_mean(),
        logger,
    );

    // Mean-Variance Modeling
    let adj_var = model_mean_variance(&filt_matrix, n_controls, config.model_choice(), logger);

    // sgRNA Ranking (Enrichment)
    let sgrna_results = enrichment_testing(&filt_matrix, &adj_var, n_controls, config.correction());

    // Build sgRNA DataFrame
    let sgrna_frame = SgrnaFrame::new(
        &filt_sgrna_names,
        &filt_gene_names,
        &adj_var,
        &sgrna_results,
    );
    sgrna_frame.write(config.prefix())?;

    if skip_agg {
        Ok(())
    } else {
        // Gene Ranking (Aggregation)
        let aggregation_results = compute_aggregation(
            config.aggregation(),
            &sgrna_results,
            &filt_gene_names,
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
}
