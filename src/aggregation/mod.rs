mod compute_aggregation;
mod rra;
mod inc;
mod utils;

pub use compute_aggregation::{GeneAggregation, compute_aggregation};
use rra::alpha_rra;
use inc::inc;
