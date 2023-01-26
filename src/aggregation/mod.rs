mod compute_aggregation;
mod mwu_inc;
mod results;
mod rra;
mod utils;

use clap::ValueEnum;
pub use compute_aggregation::{compute_aggregation, GeneAggregation};
use mwu_inc::inc;
pub use results::AggregationResult;
use rra::alpha_rra;

/// Enum describing aggregation procedure selection
#[derive(ValueEnum, Clone, Debug)]
pub enum GeneAggregationSelection {
    /// Alpha Robust Rank Algorithm (αRRA)
    RRA,

    /// INC Method, i.e. Mann-Whitney U-Test
    Inc,
}
