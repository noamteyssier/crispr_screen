use ndarray::Array2;
use crate::utils::logging::Logger;

use super::{median_ratio_normalization, total_normalization};

#[derive(Debug)]
pub enum Normalization {
    MedianRatio,
    Total
}

/// Normalize read counts using the provided method
pub fn normalize_counts(
    count_matrix: &Array2<f64>,
    normalization: &Normalization,
    logger: &Logger) -> Array2<f64>
{
    logger.norm_method(normalization);
    match normalization{
        Normalization::MedianRatio => median_ratio_normalization(count_matrix),
        Normalization::Total => total_normalization(count_matrix)
    }
}
