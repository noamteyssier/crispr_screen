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
use bon::{bon, builder};
use geopagg::{GeoPAGG, TransformConfig, WeightConfig};
use intc::{fdr::Direction, Inc};
use ndarray::Array1;

#[builder]
struct RunAggregation<'a> {
    pvalue_low: &'a Array1<f64>,
    pvalue_high: &'a Array1<f64>,
    logfc: &'a Array1<f64>,
    gene_names: &'a Vec<String>,
    seed: u64,
    logger: &'a Logger,
}
#[bon]
impl<'a> RunAggregation<'a> {
    #[builder]
    pub fn run_rra(
        &self,
        alpha: f64,
        adjust_alpha: bool,
        npermutations: usize,
        correction: Procedure,
    ) -> InternalAggregationResult {
        let (alpha_low, alpha_high) =
            set_alpha_threshold(self.pvalue_low, self.pvalue_high, alpha, adjust_alpha);
        self.logger
            .report_rra_params(alpha_low, alpha_high, self.seed as usize);

        // Calculates the RRA score for the depleted pvalues
        let alpha_rra_low = AlphaRRA::new(
            self.gene_names,
            alpha_low,
            npermutations,
            correction,
            self.seed,
        );
        let permutation_sizes_low = alpha_rra_low
            .permutation_vectors()
            .keys()
            .copied()
            .collect::<Vec<usize>>();
        self.logger.permutation_sizes(&permutation_sizes_low);
        let result_low = alpha_rra_low
            .run(self.pvalue_low)
            .expect("Error in RRA fit for depleted pvalues");

        // Calculates the RRA score for the enriched pvalues
        let alpha_rra_high = AlphaRRA::new(
            self.gene_names,
            alpha_high,
            npermutations,
            correction,
            self.seed + 1,
        );
        let permutation_sizes_high = alpha_rra_high
            .permutation_vectors()
            .keys()
            .copied()
            .collect::<Vec<usize>>();
        self.logger.permutation_sizes(&permutation_sizes_high);
        let result_high = alpha_rra_high
            .run(self.pvalue_high)
            .expect("Error in RRA fit for enriched pvalues");

        let gene_fc_hashmap = aggregate_fold_changes(self.gene_names, self.logfc);
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

    #[builder]
    pub fn run_inc(
        &self,
        token: &str,
        fdr: f64,
        group_size: usize,
        n_draws: usize,
        num_genes: usize,
        use_product: bool,
    ) -> InternalAggregationResult {
        self.logger.report_inc_params(
            token,
            num_genes,
            fdr,
            group_size,
            n_draws,
            self.seed as usize,
        );

        let (dir_low, dir_high) = if use_product {
            (Some(Direction::Less), Some(Direction::Greater))
        } else {
            (None, None)
        };

        let result_low = Inc::new(
            self.pvalue_low,
            self.logfc,
            self.gene_names,
            token,
            num_genes,
            group_size,
            n_draws,
            fdr,
            intc::mwu::Alternative::Less,
            true,
            dir_low,
            Some(self.seed),
        )
        .fit()
        .expect("Error calculating INC on low pvalues");
        self.logger
            .report_inc_low_threshold(result_low.threshold(), use_product);

        let result_high = Inc::new(
            self.pvalue_high,
            self.logfc,
            self.gene_names,
            token,
            num_genes,
            group_size,
            n_draws,
            fdr,
            intc::mwu::Alternative::Less,
            true,
            dir_high,
            Some(self.seed),
        )
        .fit()
        .expect("Error calculating INC on high pvalues");
        self.logger
            .report_inc_high_threshold(result_high.threshold(), use_product);

        self.logger.report_inc_ntc_std(result_low.null_stddev());

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

    #[builder]
    pub fn run_geopagg(
        &self,
        token: Option<&str>,
        weight_config: WeightConfig,
        fdr: f64,
    ) -> InternalAggregationResult {
        self.logger
            .report_geopagg_params(token, fdr, weight_config, self.seed as usize);

        let geo_low = GeoPAGG::new(
            self.pvalue_low.as_slice().unwrap(),
            self.logfc.as_slice().unwrap(),
            self.gene_names,
            token,
            weight_config,
            TransformConfig::Fdr,
            self.seed as usize,
        );
        let geo_high = GeoPAGG::new(
            self.pvalue_high.as_slice().unwrap(),
            self.logfc.as_slice().unwrap(),
            self.gene_names,
            token,
            weight_config,
            TransformConfig::Fdr,
            self.seed as usize,
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
}

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

    let runner = RunAggregation::builder()
        .pvalue_low(&passing_sgrna_pvalues_low)
        .pvalue_high(&passing_sgrna_pvalues_high)
        .logfc(&passing_sgrna_logfc)
        .gene_names(&passing_gene_names)
        .seed(seed)
        .logger(logger)
        .build();

    let agg_result = match agg {
        GeneAggregation::AlpaRRA {
            alpha,
            npermutations,
            adjust_alpha,
            fdr: _,
        } => runner
            .run_rra()
            .alpha(*alpha)
            .npermutations(*npermutations)
            .adjust_alpha(*adjust_alpha)
            .correction(correction)
            .call(),

        GeneAggregation::Inc {
            token,
            fdr,
            group_size,
            n_draws,
            use_product,
        } => runner
            .run_inc()
            .token(token)
            .fdr(*fdr)
            .group_size(*group_size)
            .n_draws(*n_draws)
            .num_genes(num_genes)
            .use_product(*use_product)
            .call(),

        GeneAggregation::GeoPAGG {
            token,
            weight_config,
            fdr,
        } => runner
            .run_geopagg()
            .maybe_token(*token)
            .fdr(*fdr)
            .weight_config(*weight_config)
            .call(),
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
