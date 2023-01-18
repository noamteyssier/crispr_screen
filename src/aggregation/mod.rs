mod compute_aggregation;
mod rra;
mod mwu_inc;
mod utils;
mod results;

use clap::ValueEnum;
pub use compute_aggregation::{GeneAggregation, compute_aggregation};
pub use results::AggregationResult;
use rra::alpha_rra;
use mwu_inc::inc;


/// Enum describing aggregation procedure selection
#[derive(ValueEnum, Clone, Debug)]
pub enum GeneAggregationSelection {

    /// Alpha Robust Rank Algorithm (αRRA)
    RRA,

    /// INC Method, i.e. Mann-Whitney U-Test
    Inc,
}
