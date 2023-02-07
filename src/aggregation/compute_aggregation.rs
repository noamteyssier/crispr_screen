use super::{
    utils::{filter_zeros, set_alpha_threshold},
    AggregationResult, GeneAggregation,
};
use crate::{
    enrich::EnrichmentResult,
    utils::{agg::aggregate_fold_changes, logging::Logger},
};
use adjustp::Procedure;
use alpha_rra::AlphaRRA;
use intc::Inc;
use ndarray::Array1;

/// Aggregates the results of the gene aggregation analysis for internal use
struct InternalAggregationResult {
    genes: Vec<String>,
    scores_low: Array1<f64>,
    pvalues_low: Array1<f64>,
    correction_low: Array1<f64>,
    scores_high: Array1<f64>,
    pvalues_high: Array1<f64>,
    correction_high: Array1<f64>,
}
impl InternalAggregationResult {
    pub fn new(
        genes: Vec<String>,
        scores_low: Array1<f64>,
        pvalues_low: Array1<f64>,
        correction_low: Array1<f64>,
        scores_high: Array1<f64>,
        pvalues_high: Array1<f64>,
        correction_high: Array1<f64>,
    ) -> Self {
        Self {
            genes,
            scores_low,
            pvalues_low,
            correction_low,
            scores_high,
            pvalues_high,
            correction_high,
        }
    }
}

/// Runs the RRA gene aggregation procedure
fn run_rra(
    pvalue_low: &Array1<f64>,
    pvalue_high: &Array1<f64>,
    gene_names: &[String],
    alpha: f64,
    adjust_alpha: bool,
    npermutations: usize,
    correction: Procedure,
    logger: &Logger,
) -> InternalAggregationResult {
    let (alpha_low, alpha_high) =
        set_alpha_threshold(&pvalue_low, &pvalue_high, alpha, adjust_alpha);
    logger.report_rra_alpha(alpha_low, alpha_high);

    // Calculates the RRA score for the depleted pvalues
    let alpha_rra_low = AlphaRRA::new(gene_names, alpha_low, npermutations, correction);
    let permutation_sizes_low = alpha_rra_low
        .permutation_vectors()
        .keys()
        .cloned()
        .collect::<Vec<usize>>();
    logger.permutation_sizes(&permutation_sizes_low);
    let result_low = alpha_rra_low
        .run(pvalue_low)
        .expect("Error in RRA fit for depleted pvalues");

    // Calculates the RRA score for the enriched pvalues
    let alpha_rra_high = AlphaRRA::new(gene_names, alpha_high, npermutations, correction);
    let permutation_sizes_high = alpha_rra_high
        .permutation_vectors()
        .keys()
        .cloned()
        .collect::<Vec<usize>>();
    logger.permutation_sizes(&permutation_sizes_high);
    let result_high = alpha_rra_high
        .run(pvalue_high)
        .expect("Error in RRA fit for enriched pvalues");

    InternalAggregationResult::new(
        result_low.names().to_vec(),
        result_low.scores().to_owned(),
        result_low.pvalues().to_owned(),
        result_low.adj_pvalues().to_owned(),
        result_high.scores().to_owned(),
        result_high.pvalues().to_owned(),
        result_high.adj_pvalues().to_owned(),
    )
}

/// Runs the INC gene aggregation procedure
fn run_inc(
    pvalue_low: &Array1<f64>,
    pvalue_high: &Array1<f64>,
    gene_names: &[String],
    token: &str,
    fdr: f64,
    group_size: usize,
    num_genes: usize,
    logger: &Logger,
) -> InternalAggregationResult {
    logger.report_inc_params(token, num_genes, fdr, group_size);
    let result_low = Inc::new(
        pvalue_low,
        gene_names,
        token,
        num_genes,
        group_size,
        fdr,
        intc::mwu::Alternative::Less,
        true,
    )
    .fit()
    .expect("Error calculating INC on low pvalues");
    logger.report_inc_low_threshold(result_low.threshold());

    let result_high = Inc::new(
        pvalue_high,
        gene_names,
        token,
        num_genes,
        group_size,
        fdr,
        intc::mwu::Alternative::Less,
        true,
    )
    .fit()
    .expect("Error calculating INC on high pvalues");
    logger.report_inc_high_threshold(result_high.threshold());

    InternalAggregationResult::new(
        result_low.genes().to_vec(),
        result_low.u_scores().to_owned(),
        result_low.u_pvalues().to_owned(),
        result_low.fdr().to_owned(),
        result_high.u_scores().to_owned(),
        result_high.u_pvalues().to_owned(),
        result_high.fdr().to_owned(),
    )
}

/// Computes gene aggregation using the provided method and associated configurations.
pub fn compute_aggregation(
    agg: &GeneAggregation,
    sgrna_results: &EnrichmentResult,
    gene_names: &[String],
    logger: &Logger,
    correction: Procedure,
) -> AggregationResult {
    logger.start_gene_aggregation();

    let gene_fc_hashmap = aggregate_fold_changes(
        gene_names,
        sgrna_results.fold_change(),
        sgrna_results.pvalues_twosided(),
    );
    let num_genes = gene_fc_hashmap.len();

    let (passing_gene_names, passing_sgrna_pvalues_low, passing_sgrna_pvalues_high) = filter_zeros(
        sgrna_results.base_means(),
        gene_names,
        sgrna_results.pvalues_low(),
        sgrna_results.pvalues_high(),
        logger,
    );

    let agg_result = match agg {
        GeneAggregation::AlpaRRA {
            alpha,
            npermutations,
            adjust_alpha,
        } => run_rra(
            &passing_sgrna_pvalues_low,
            &passing_sgrna_pvalues_high,
            &passing_gene_names,
            *alpha,
            *adjust_alpha,
            *npermutations,
            correction,
            logger,
        ),
        GeneAggregation::Inc {
            token,
            fdr,
            group_size,
        } => run_inc(
            &passing_sgrna_pvalues_low,
            &passing_sgrna_pvalues_high,
            &passing_gene_names,
            token,
            *fdr,
            *group_size,
            num_genes,
            logger,
        ),
    };

    let gene_fc = agg_result
        .genes
        .iter()
        .map(|x| gene_fc_hashmap.get(x).unwrap_or(&0.0))
        .copied()
        .collect::<Array1<f64>>();

    AggregationResult::new(
        agg_result.genes,
        gene_fc,
        agg_result.pvalues_low,
        agg_result.pvalues_high,
        agg_result.correction_low,
        agg_result.correction_high,
        agg_result.scores_low,
        agg_result.scores_high,
    )
}
