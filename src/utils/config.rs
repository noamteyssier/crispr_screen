use crate::{
    aggregation::GeneAggregation, enrich::TestStrategy, model::ModelChoice, norm::Normalization,
};
use adjustp::Procedure;
use getset::Getters;

#[derive(Debug, Getters)]
#[getset(get = "pub")]
pub struct Configuration<'a> {
    normalization: Normalization,
    aggregation: GeneAggregation<'a>,
    correction: Procedure,
    model_choice: ModelChoice,
    min_base_mean: f64,
    strategy: TestStrategy,
    seed: u64,
    prefix: &'a str,
}
impl<'a> Configuration<'a> {
    pub fn new(
        normalization: Normalization,
        aggregation: GeneAggregation<'a>,
        correction: Procedure,
        model_choice: ModelChoice,
        min_base_mean: f64,
        strategy: TestStrategy,
        seed: u64,
        prefix: &'a str,
    ) -> Self {
        Self {
            normalization,
            aggregation,
            correction,
            model_choice,
            min_base_mean,
            strategy,
            seed,
            prefix,
        }
    }
    pub fn new_agg(
        aggregation: GeneAggregation<'a>,
        correction: Procedure,
        seed: u64,
        prefix: &'a str,
    ) -> Self {
        Self {
            normalization: Normalization::default(),
            aggregation,
            correction,
            model_choice: ModelChoice::default(),
            min_base_mean: 0.0,
            strategy: TestStrategy::default(),
            seed,
            prefix,
        }
    }
}

#[cfg(test)]
mod testing {
    use adjustp::Procedure;

    use crate::{aggregation::GeneAggregation, model::ModelChoice, norm::Normalization};

    use super::Configuration;

    fn build_config<'a>() -> Configuration<'a> {
        let normalization = Normalization::MedianRatio;
        let aggregation = GeneAggregation::Inc {
            token: "non-targeting",
            fdr: 0.05,
            group_size: 5,
            use_product: true,
            n_draws: 100,
        };
        let correction = Procedure::BenjaminiHochberg;
        let model_choice = ModelChoice::Wols;
        let base_mean = 0.0;
        let seed = 0;
        let prefix = "results";
        Configuration::new(
            normalization,
            aggregation,
            correction,
            model_choice,
            base_mean,
            seed,
            prefix,
        )
    }

    #[test]
    fn test_normalization() {
        let config = build_config();
        assert_eq!(config.normalization(), &Normalization::MedianRatio);
    }

    #[test]
    fn test_aggregation() {
        let config = build_config();
        match config.aggregation() {
            GeneAggregation::Inc {
                token,
                fdr,
                group_size,
                use_product,
                n_draws,
            } => {
                assert_eq!(token, &"non-targeting");
                assert_eq!(fdr, &0.05);
                assert_eq!(group_size, &5);
                assert_eq!(use_product, &true);
                assert_eq!(n_draws, &100);
            }
            _ => panic!(),
        }
    }

    #[test]
    fn test_correction() {
        let config = build_config();
        assert_eq!(config.correction(), Procedure::BenjaminiHochberg);
    }

    #[test]
    fn test_model_choice() {
        let config = build_config();
        assert_eq!(*config.model_choice(), ModelChoice::Wols);
    }

    #[test]
    fn test_prefix() {
        let config = build_config();
        assert_eq!(config.prefix(), "results");
    }
}
