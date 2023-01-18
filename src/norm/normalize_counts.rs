use clap::ValueEnum;
use ndarray::Array2;
use super::{median_ratio_normalization, total_normalization};

#[derive(ValueEnum, Debug, Clone)]
pub enum Normalization {
    MedianRatio,
    Total
}

/// Normalize read counts using the provided method
pub fn normalize_counts(
    count_matrix: &Array2<f64>,
    normalization: &Normalization) -> Array2<f64>
{
    match normalization{
        Normalization::MedianRatio => median_ratio_normalization(count_matrix),
        Normalization::Total => total_normalization(count_matrix)
    }
}
