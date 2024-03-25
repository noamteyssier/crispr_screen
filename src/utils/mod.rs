pub mod agg;
pub mod config;
pub mod filter;
pub mod logging;
pub mod math;
use clap::ValueEnum;

/// Selection for P-Value Adjustments
#[derive(ValueEnum, Clone, Debug)]
pub enum Adjustment {
    /// Bonferroni
    Bf,

    /// Benjamini-Hochberg
    Bh,

    /// Benjamini-Yekutieli
    By,
}
