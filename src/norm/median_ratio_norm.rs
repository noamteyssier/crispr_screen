use anyhow::{bail, Result};
use ndarray::{Array1, Array2, ArrayView1, Axis};
use ndarray_stats::SummaryStatisticsExt;
use std::ops::{Div, Mul};

/// Calculates the median of a provided ndarray
pub fn median(array: &ArrayView1<f64>) -> f64 {
    let mut sorted = array.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("NaN Uncovered in Median"));
    if array.len() % 2 == 0 {
        let rhs = array.len().div(2);
        let lhs = rhs - 1;
        (sorted[lhs] + sorted[rhs]).div(2.)
    } else {
        let midpoint = array.len().div(2);
        sorted[midpoint]
    }
}

/// Calculates the geometric means of each row in a 2D matrix
/// This will also set `0.`s to `1.`s to avoid downstream numerical instability
fn geometric_means(matrix: &Array2<f64>) -> Array1<f64> {
    matrix
        .map_axis(Axis(1), |row| {
            row.geometric_mean().expect("Unexpected Empty Row")
        })
        .iter()
        .map(|x| if *x == 0. { 1. } else { *x })
        .collect()
}

/// Performs the median ratio normalization method.
/// Read counts are adjusted by the median of the size factors of each `sgRNA`.
/// Size factors are computed as the geometric mean of the `sgRNAs` across all
/// experimental libraries.
pub fn median_ratio_normalization(matrix: &Array2<f64>) -> Result<Array2<f64>> {
    let gmeans = geometric_means(matrix);
    let transformed_matrix = matrix.t().div(gmeans).reversed_axes();
    let sample_scalars = transformed_matrix.map_axis(Axis(0), |axis| 1. / median(&axis));

    if sample_scalars.iter().any(|x| x.is_infinite()) {
        bail!("median-ratio is unstable")
    }

    Ok(matrix.mul(sample_scalars))
}

#[cfg(test)]
mod testing {
    use super::{median, median_ratio_normalization};
    use ndarray::{Array1, Array2};
    use ndarray_rand::{rand_distr::Uniform, RandomExt};

    #[test]
    fn test_median_even() {
        let arr = Array1::range(1., 5., 1.);
        assert_eq!(median(&arr.view()), 2.5);
    }

    #[test]
    fn test_median_odd() {
        let arr = Array1::range(1., 6., 1.);
        assert_eq!(median(&arr.view()), 3.0);
    }

    #[test]
    fn test_median_normalization() {
        (0..1000).for_each(|_| {
            let matrix = Array2::random((10, 4), Uniform::new(1., 10.));
            let norm = median_ratio_normalization(&matrix).unwrap();

            // matrices must be equal shape
            assert_eq!(norm.shape(), matrix.shape());
        })
    }

    #[test]
    fn test_median_error() {
        let matrix = Array2::zeros((10, 4));
        let norm = median_ratio_normalization(&matrix);
        assert!(norm.is_err());
    }
}
