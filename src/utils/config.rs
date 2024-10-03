use crate::{
    aggregation::GeneAggregation, enrich::TestStrategy, model::ModelChoice, norm::Normalization,
};
use adjustp::Procedure;
use bon::{builder, Builder};
use getset::Getters;

#[derive(Debug, Getters, Builder)]
#[getset(get = "pub")]
pub struct Configuration<'a> {
    #[builder(default)]
    normalization: Normalization,
    aggregation: GeneAggregation<'a>,
    #[builder(default = Procedure::BenjaminiHochberg)]
    correction: Procedure,
    #[builder(default)]
    model_choice: ModelChoice,
    #[builder(default)]
    min_base_mean: f64,
    #[builder(default)]
    strategy: TestStrategy,
    #[builder(default)]
    seed: u64,
    prefix: &'a str,
}

#[cfg(test)]
mod testing {
    use adjustp::Procedure;

    use crate::{
        aggregation::GeneAggregation, enrich::TestStrategy, model::ModelChoice, norm::Normalization,
    };

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
        let strategy = TestStrategy::default();
        let seed = 0;
        let prefix = "results";
        Configuration::builder()
            .normalization(normalization)
            .aggregation(aggregation)
            .correction(correction)
            .model_choice(model_choice)
            .min_base_mean(base_mean)
            .strategy(strategy)
            .seed(seed)
            .prefix(prefix)
            .build()
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
        assert_eq!(*config.correction(), Procedure::BenjaminiHochberg);
    }

    #[test]
    fn test_model_choice() {
        let config = build_config();
        assert_eq!(*config.model_choice(), ModelChoice::Wols);
    }

    #[test]
    fn test_prefix() {
        let config = build_config();
        assert_eq!(*config.prefix(), "results");
    }
}
