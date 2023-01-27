mod compute_aggregation;
mod mwu_inc;
mod results;
mod rra;
mod utils;

use clap::ValueEnum;
pub use compute_aggregation::compute_aggregation;
use mwu_inc::inc;
pub use results::AggregationResult;
use rra::alpha_rra;

/// Enum describing aggregation procedure selection
#[derive(ValueEnum, Clone, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub enum GeneAggregationSelection {
    /// Alpha Robust Rank Algorithm (αRRA)
    RRA,

    /// INC Method, i.e. Mann-Whitney U-Test
    Inc,
}

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
