use crate::utils::logging::Logger;
use clap::ValueEnum;
use ndarray::Array2;

use super::{median_ratio_normalization, total_normalization};

#[derive(ValueEnum, Debug, Clone)]
pub enum Normalization {
    /// Median ratio of sgRNA geometric means (unstable at high zeros)
    MedianRatio,

    /// Total read count per sample scaling (more stable)
    Total,
}

/// Normalize read counts using the provided method
pub fn normalize_counts(
    count_matrix: &Array2<f64>,
    normalization: &Normalization,
    logger: &Logger,
) -> Array2<f64> {
    match normalization {
        Normalization::MedianRatio => {
            if let Ok(normed_matrix) = median_ratio_normalization(count_matrix) {
                normed_matrix
            } else {
                logger.convert_normalization();
                total_normalization(count_matrix)
            }
        }
        Normalization::Total => total_normalization(count_matrix),
    }
}
