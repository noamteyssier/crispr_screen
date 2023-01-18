mod dfutils;
pub mod io;
pub mod math;
pub mod logging;
use clap::ValueEnum;
pub use dfutils::{parse_to_string_vec, parse_to_ndarray, parse_genes, parse_sgrna, FormatError};

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
