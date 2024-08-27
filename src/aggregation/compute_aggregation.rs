use super::{
    utils::{filter_zeros, num_unique, set_alpha_threshold},
    AggregationResult, GeneAggregation,
};
use crate::{
    enrich::EnrichmentResult,
    utils::{agg::aggregate_fold_changes, logging::Logger},
};
use adjustp::Procedure;
use alpha_rra::AlphaRRA;
use geopagg::{GeoPAGG, TransformConfig, WeightConfig};
use intc::{fdr::Direction, Inc};
use ndarray::Array1;

/// Aggregates the results of the gene aggregation analysis for internal use
struct InternalAggregationResult {
    genes: Vec<String>,
    logfc: Array1<f64>,
    scores_low: Array1<f64>,
    pvalues_low: Array1<f64>,
    correction_low: Array1<f64>,
    scores_high: Array1<f64>,
    pvalues_high: Array1<f64>,
    correction_high: Array1<f64>,
    threshold_low: Option<f64>,
    threshold_high: Option<f64>,
}
impl InternalAggregationResult {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        genes: Vec<String>,
        logfc: Array1<f64>,
        scores_low: Array1<f64>,
        pvalues_low: Array1<f64>,
        correction_low: Array1<f64>,
        scores_high: Array1<f64>,
        pvalues_high: Array1<f64>,
        correction_high: Array1<f64>,
        threshold_low: Option<f64>,
        threshold_high: Option<f64>,
    ) -> Self {
        Self {
            genes,
            logfc,
            scores_low,
            pvalues_low,
            correction_low,
            scores_high,
            pvalues_high,
            correction_high,
            threshold_low,
            threshold_high,
        }
    }
}

/// Runs the RRA gene aggregation procedure
#[allow(clippy::too_many_arguments)]
fn run_rra(
    pvalue_low: &Array1<f64>,
    pvalue_high: &Array1<f64>,
    logfc: &Array1<f64>,
    gene_names: &[String],
    alpha: f64,
    adjust_alpha: bool,
    npermutations: usize,
    correction: Procedure,
    seed: u64,
    logger: &Logger,
) -> InternalAggregationResult {
    let (alpha_low, alpha_high) = set_alpha_threshold(pvalue_low, pvalue_high, alpha, adjust_alpha);
    logger.report_rra_params(alpha_low, alpha_high, seed as usize);

    // Calculates the RRA score for the depleted pvalues
    let alpha_rra_low = AlphaRRA::new(gene_names, alpha_low, npermutations, correction, seed);
    let permutation_sizes_low = alpha_rra_low
        .permutation_vectors()
        .keys()
        .copied()
        .collect::<Vec<usize>>();
    logger.permutation_sizes(&permutation_sizes_low);
    let result_low = alpha_rra_low
        .run(pvalue_low)
        .expect("Error in RRA fit for depleted pvalues");

    // Calculates the RRA score for the enriched pvalues
    let alpha_rra_high = AlphaRRA::new(gene_names, alpha_high, npermutations, correction, seed + 1);
    let permutation_sizes_high = alpha_rra_high
        .permutation_vectors()
        .keys()
        .copied()
        .collect::<Vec<usize>>();
    logger.permutation_sizes(&permutation_sizes_high);
    let result_high = alpha_rra_high
        .run(pvalue_high)
        .expect("Error in RRA fit for enriched pvalues");

    let gene_fc_hashmap = aggregate_fold_changes(gene_names, logfc);
    let gene_fc = result_low
        .names()
        .iter()
        .map(|gene| gene_fc_hashmap.get(gene).unwrap_or(&0.0))
        .copied()
        .collect();

    InternalAggregationResult::new(
        result_low.names().to_vec(),
        gene_fc,
        result_low.scores().to_owned(),
        result_low.pvalues().to_owned(),
        result_low.adj_pvalues().to_owned(),
        result_high.scores().to_owned(),
        result_high.pvalues().to_owned(),
        result_high.adj_pvalues().to_owned(),
        None,
        None,
    )
}

/// Runs the INC gene aggregation procedure
#[allow(clippy::too_many_arguments)]
fn run_inc(
    pvalue_low: &Array1<f64>,
    pvalue_high: &Array1<f64>,
    log2_fold_change: &Array1<f64>,
    gene_names: &[String],
    token: &str,
    fdr: f64,
    group_size: usize,
    n_draws: usize,
    num_genes: usize,
    use_product: bool,
    seed: u64,
    logger: &Logger,
) -> InternalAggregationResult {
    logger.report_inc_params(token, num_genes, fdr, group_size, n_draws, seed as usize);

    let (dir_low, dir_high) = if use_product {
        (Some(Direction::Less), Some(Direction::Greater))
    } else {
        (None, None)
    };

    let result_low = Inc::new(
        pvalue_low,
        log2_fold_change,
        gene_names,
        token,
        num_genes,
        group_size,
        n_draws,
        fdr,
        intc::mwu::Alternative::Less,
        true,
        dir_low,
        Some(seed),
    )
    .fit()
    .expect("Error calculating INC on low pvalues");
    logger.report_inc_low_threshold(result_low.threshold(), use_product);

    let result_high = Inc::new(
        pvalue_high,
        log2_fold_change,
        gene_names,
        token,
        num_genes,
        group_size,
        n_draws,
        fdr,
        intc::mwu::Alternative::Less,
        true,
        dir_high,
        Some(seed),
    )
    .fit()
    .expect("Error calculating INC on high pvalues");
    logger.report_inc_high_threshold(result_high.threshold(), use_product);

    logger.report_inc_ntc_std(result_low.null_stddev());

    InternalAggregationResult::new(
        result_low.genes().to_vec(),
        result_low.logfc().to_owned(),
        result_low.u_scores().to_owned(),
        result_low.u_pvalues().to_owned(),
        result_low.fdr().to_owned(),
        result_high.u_scores().to_owned(),
        result_high.u_pvalues().to_owned(),
        result_high.fdr().to_owned(),
        Some(result_low.threshold()),
        Some(result_high.threshold()),
    )
}

fn run_geopagg(
    pvalue_low: &Array1<f64>,
    pvalue_high: &Array1<f64>,
    logfc: &Array1<f64>,
    gene_names: &[String],
    token: Option<&str>,
    weight_config: WeightConfig,
    fdr: f64,
    seed: u64,
    logger: &Logger,
) -> InternalAggregationResult {
    logger.report_geopagg_params(token, fdr, weight_config, seed as usize);

    let geo_low = GeoPAGG::new(
        pvalue_low.as_slice().unwrap(),
        logfc.as_slice().unwrap(),
        gene_names,
        token,
        weight_config,
        TransformConfig::Fdr,
        seed as usize,
    );
    let geo_high = GeoPAGG::new(
        pvalue_high.as_slice().unwrap(),
        logfc.as_slice().unwrap(),
        gene_names,
        token,
        weight_config,
        TransformConfig::Fdr,
        seed as usize,
    );

    let geo_low_results = geo_low.run();
    let geo_high_results = geo_high.run();

    InternalAggregationResult::new(
        geo_low_results.genes.to_vec(),
        Array1::from(geo_low_results.logfcs),
        Array1::from(geo_low_results.empirical_fdr),
        Array1::from(geo_low_results.wgms),
        Array1::from(geo_low_results.adjusted_empirical_fdr),
        Array1::from(geo_high_results.empirical_fdr),
        Array1::from(geo_high_results.wgms),
        Array1::from(geo_high_results.adjusted_empirical_fdr),
        Some(fdr),
        Some(fdr),
    )
}

/// Computes gene aggregation using the provided method and associated configurations.
pub fn compute_aggregation(
    agg: &GeneAggregation,
    sgrna_results: &EnrichmentResult,
    gene_names: &[String],
    logger: &Logger,
    correction: Procedure,
    seed: u64,
) -> AggregationResult {
    logger.start_gene_aggregation();

    let num_genes = num_unique(gene_names);

    let (
        passing_gene_names,
        passing_sgrna_pvalues_low,
        passing_sgrna_pvalues_high,
        passing_sgrna_logfc,
    ) = filter_zeros(
        sgrna_results.base_means(),
        gene_names,
        sgrna_results.pvalues_low(),
        sgrna_results.pvalues_high(),
        sgrna_results.log_fold_change(),
        logger,
    );

    let agg_result = match agg {
        GeneAggregation::AlpaRRA {
            alpha,
            npermutations,
            adjust_alpha,
            fdr: _,
        } => run_rra(
            &passing_sgrna_pvalues_low,
            &passing_sgrna_pvalues_high,
            &passing_sgrna_logfc,
            &passing_gene_names,
            *alpha,
            *adjust_alpha,
            *npermutations,
            correction,
            seed,
            logger,
        ),
        GeneAggregation::Inc {
            token,
            fdr,
            group_size,
            n_draws,
            use_product,
        } => run_inc(
            &passing_sgrna_pvalues_low,
            &passing_sgrna_pvalues_high,
            &passing_sgrna_logfc,
            &passing_gene_names,
            token,
            *fdr,
            *group_size,
            *n_draws,
            num_genes,
            *use_product,
            seed,
            logger,
        ),
        GeneAggregation::GeoPAGG {
            token,
            weight_config,
            fdr,
        } => run_geopagg(
            &passing_sgrna_pvalues_low,
            &passing_sgrna_pvalues_high,
            &passing_sgrna_logfc,
            &passing_gene_names,
            *token,
            *weight_config,
            *fdr,
            seed,
            logger,
        ),
    };

    let fold_change = agg_result
        .logfc
        .iter()
        .map(|x| x.exp2())
        .collect::<Array1<f64>>();

    AggregationResult::new(
        agg_result.genes,
        fold_change,
        agg_result.pvalues_low,
        agg_result.pvalues_high,
        agg_result.correction_low,
        agg_result.correction_high,
        agg_result.scores_low,
        agg_result.scores_high,
        agg_result.threshold_low,
        agg_result.threshold_high,
    )
}
