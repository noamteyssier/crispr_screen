pub mod math;
pub mod logging;
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
