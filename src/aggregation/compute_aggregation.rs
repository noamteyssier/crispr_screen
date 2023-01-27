use super::{alpha_rra, inc, AggregationResult, utils::{set_alpha_threshold, filter_zeros}};
use crate::{
    enrich::EnrichmentResult,
    utils::{agg::aggregate_fold_changes, logging::Logger},
};
use adjustp::Procedure;
use ndarray::{Array1, Array2};

/// Enum describing the different gene aggregation procedures and their associated configurations.
#[derive(Debug)]
pub enum GeneAggregation<'a> {
    AlpaRRA {
        alpha: f64,
        npermutations: usize,
        adjust_alpha: bool,
    },
    Inc {
        token: &'a str,
    },
}

/// Aggregates the results of the gene aggregation analysis for internal use
struct InternalAggregationResult {
    genes: Vec<String>,
    gene_scores_low: Array1<f64>,
    gene_pvalues_low: Array1<f64>,
    gene_scores_high: Array1<f64>,
    gene_pvalues_high: Array1<f64>,
}
impl InternalAggregationResult {
    pub fn new(
        genes: Vec<String>,
        gene_scores_low: Array1<f64>,
        gene_pvalues_low: Array1<f64>,
        gene_scores_high: Array1<f64>,
        gene_pvalues_high: Array1<f64>,
    ) -> Self {
        Self {
            genes,
            gene_scores_low,
            gene_pvalues_low,
            gene_scores_high,
            gene_pvalues_high,
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
    logger: &Logger,
) -> InternalAggregationResult {
    let (alpha_low, alpha_high) = set_alpha_threshold(
        &pvalue_low,
        &pvalue_high,
        alpha,
        adjust_alpha,
    );
    logger.report_rra_alpha(alpha_low, alpha_high);
    let (genes, gene_scores_low, gene_pvalues_low) = alpha_rra(
        pvalue_low,
        gene_names,
        alpha_low,
        npermutations,
        logger,
    );
    let (_, gene_scores_high, gene_pvalues_high) = alpha_rra(
        pvalue_high,
        gene_names,
        alpha_high,
        npermutations,
        logger,
    );
    InternalAggregationResult::new(genes, gene_scores_low, gene_pvalues_low, gene_scores_high, gene_pvalues_high)
}

/// Runs the INC gene aggregation procedure
fn run_inc(
    pvalue_low: &Array1<f64>,
    pvalue_high: &Array1<f64>,
    gene_names: &[String],
    token: &str,
    logger: &Logger,
) -> InternalAggregationResult {
    let (genes, gene_scores_low, gene_pvalues_low) = inc(
        pvalue_low,
        gene_names,
        token,
        logger,
    );
    let (_, gene_scores_high, gene_pvalues_high) = inc(
        pvalue_high,
        gene_names,
        token,
        logger,
    );
    InternalAggregationResult::new(genes, gene_scores_low, gene_pvalues_low, gene_scores_high, gene_pvalues_high)
}


/// Computes gene aggregation using the provided method and associated configurations.
pub fn compute_aggregation(
    agg: &GeneAggregation,
    normed_matrix: &Array2<f64>,
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

    let (passing_gene_names, passing_sgrna_pvalues_low, passing_sgrna_pvalues_high) = filter_zeros(
        normed_matrix,
        gene_names,
        sgrna_results.pvalues_low(),
        sgrna_results.pvalues_high(),
        logger,
    );

    let agg_result = match agg
    {
        GeneAggregation::AlpaRRA {
            alpha,
            npermutations,
            adjust_alpha,
        } => {
            run_rra(
                &passing_sgrna_pvalues_low,
                &passing_sgrna_pvalues_high,
                &passing_gene_names,
                *alpha,
                *adjust_alpha,
                *npermutations,
                logger,
            )

        }
        GeneAggregation::Inc { token } => {
            run_inc(
                &passing_sgrna_pvalues_low,
                &passing_sgrna_pvalues_high,
                &passing_gene_names,
                token,
                logger,
            )
        }
    };

    let gene_fc = agg_result.genes
        .iter()
        .map(|x| gene_fc_hashmap.get(x).unwrap_or(&0.0))
        .copied()
        .collect::<Array1<f64>>();

    AggregationResult::new(
        agg_result.genes,
        gene_fc,
        agg_result.gene_pvalues_low,
        agg_result.gene_pvalues_high,
        agg_result.gene_scores_low,
        agg_result.gene_scores_high,
        correction,
    )
}
