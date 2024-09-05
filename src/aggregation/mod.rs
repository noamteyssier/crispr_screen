mod compute_aggregation;
mod results;
mod utils;

use clap::ValueEnum;
pub use compute_aggregation::compute_aggregation;
use geopagg::WeightConfig;
pub use results::AggregationResult;

/// Enum describing aggregation procedure selection
#[derive(ValueEnum, Clone, Debug, PartialEq)]
#[allow(clippy::upper_case_acronyms)]
pub enum GeneAggregationSelection {
    /// Alpha Robust Rank Algorithm (αRRA)
    RRA,

    /// INC Method, i.e. Mann-Whitney U-Test
    Inc,

    /// GeoPAGG Method
    #[value(name = "geopagg")]
    GeoPAGG,
}

/// Enum describing the different gene aggregation procedures and their associated configurations.
#[derive(Debug)]
pub enum GeneAggregation<'a> {
    AlpaRRA {
        alpha: f64,
        npermutations: usize,
        adjust_alpha: bool,
        fdr: f64,
    },
    Inc {
        token: &'a str,
        fdr: f64,
        group_size: usize,
        n_draws: usize,
        use_product: bool,
    },
    GeoPAGG {
        token: Option<&'a str>,
        weight_config: WeightConfig,
        fdr: f64,
        use_product: bool,
    },
}

#[derive(ValueEnum, Clone, Copy, Debug, Default)]
pub enum GeoPAGGWeightConfigEnum {
    #[default]
    DropFirst,
    RankOrder,
    Balanced,
}

#[cfg(test)]
mod testing {
    use super::{GeneAggregation, GeneAggregationSelection};
    use clap::ValueEnum;

    #[test]
    fn test_enum() {
        let rra = GeneAggregationSelection::RRA;
        let inc = GeneAggregationSelection::Inc;

        assert_eq!(
            GeneAggregationSelection::from_str("RRA", true).unwrap(),
            rra
        );
        assert_eq!(
            GeneAggregationSelection::from_str("Inc", true).unwrap(),
            inc
        );
    }

    #[test]
    fn test_enum_case_insensitive() {
        let rra = GeneAggregationSelection::RRA;
        let inc = GeneAggregationSelection::Inc;

        assert_eq!(
            GeneAggregationSelection::from_str("rra", true).unwrap(),
            rra
        );
        assert_eq!(
            GeneAggregationSelection::from_str("inc", true).unwrap(),
            inc
        );
    }

    #[test]
    fn test_enum_invalid() {
        assert!(GeneAggregationSelection::from_str("invalid", true).is_err());
    }

    #[test]
    fn test_enum_rra_string() {
        let rra = GeneAggregationSelection::RRA;
        assert_eq!(format!("{:?}", rra), "RRA");
    }

    #[test]
    fn test_enum_inc_string() {
        let inc = GeneAggregationSelection::Inc;
        assert_eq!(format!("{:?}", inc), "Inc");
    }

    #[test]
    fn test_enum_gene_aggregation_rra() {
        let agg = GeneAggregation::AlpaRRA {
            alpha: 0.5,
            npermutations: 100,
            adjust_alpha: true,
            fdr: 0.05,
        };
        assert_eq!(
            format!("{:?}", agg),
            "AlpaRRA { alpha: 0.5, npermutations: 100, adjust_alpha: true, fdr: 0.05 }"
        );
    }

    #[test]
    fn test_enum_gene_aggregation_inc() {
        let agg = GeneAggregation::Inc {
            token: "non-targeting",
            fdr: 0.05,
            group_size: 3,
            n_draws: 100,
            use_product: true,
        };
        assert_eq!(
            format!("{:?}", agg),
            "Inc { token: \"non-targeting\", fdr: 0.05, group_size: 3, n_draws: 100, use_product: true }"
        );
    }
}
