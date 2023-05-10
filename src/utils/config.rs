use crate::{aggregation::GeneAggregation, model::ModelChoice, norm::Normalization};
use adjustp::Procedure;

#[derive(Debug)]
pub struct Configuration<'a> {
    normalization: Normalization,
    aggregation: GeneAggregation<'a>,
    correction: Procedure,
    model_choice: ModelChoice,
    seed: u64,
    prefix: &'a str,
}
impl<'a> Configuration<'a> {
    pub fn new(
        normalization: Normalization,
        aggregation: GeneAggregation<'a>,
        correction: Procedure,
        model_choice: ModelChoice,
        seed: u64,
        prefix: &'a str,
    ) -> Self {
        Self {
            normalization,
            aggregation,
            correction,
            model_choice,
            seed,
            prefix,
        }
    }
    pub fn normalization(&self) -> &Normalization {
        &self.normalization
    }
    pub fn aggregation(&self) -> &GeneAggregation {
        &self.aggregation
    }
    pub fn correction(&self) -> Procedure {
        self.correction
    }
    pub fn model_choice(&self) -> &ModelChoice {
        &self.model_choice
    }
    pub fn seed(&self) -> u64 {
        self.seed
    }
    pub fn prefix(&self) -> &str {
        self.prefix
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
        };
        let correction = Procedure::BenjaminiHochberg;
        let model_choice = ModelChoice::Wols;
        let seed = 0;
        let prefix = "results";
        Configuration::new(normalization, aggregation, correction, model_choice, seed, prefix)
    }

    #[test]
    fn test_normalization() {
        let config = build_config();
        match config.normalization() {
            Normalization::MedianRatio => assert!(true),
            _ => assert!(false),
        }
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
            } => {
                assert_eq!(token, &"non-targeting");
                assert_eq!(fdr, &0.05);
                assert_eq!(group_size, &5);
            }
            _ => assert!(false),
        }
    }

    #[test]
    fn test_correction() {
        let config = build_config();
        match config.correction {
            Procedure::BenjaminiHochberg => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    fn test_model_choice() {
        let config = build_config();
        match config.model_choice() {
            ModelChoice::Wols => assert!(true),
            _ => assert!(false),
        }
    }

    #[test]
    fn test_prefix() {
        let config = build_config();
        assert_eq!(config.prefix(), "results");
    }
}
