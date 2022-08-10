use std::ops::{Div, Mul};
use ndarray::{Axis, Array2, ArrayView1, Array1};
use ndarray_stats::SummaryStatisticsExt;

/// Calculates the median of a provided ndarray
fn median(
    array: &ArrayView1<f64>) -> f64
{
    let mut sorted = array.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("NaN Uncovered in Median"));
    if array.len() % 2 == 0 {
        let midpoint = array.len().div(2);
        sorted[midpoint]
    } else {
        let lhs = array.len().div(2);
        let rhs = lhs + 1;
        (sorted[lhs] + sorted[rhs]).div(2.)
    }
}

/// Calculates the geometric means of each row in a 2D matrix
/// This will also set `0.`s to `1.`s to avoid downstream numerical instability 
fn geometric_means(matrix: &Array2<f64>) -> Array1<f64>
{
    matrix
        .map_axis(
            Axis(1), |row| row.geometric_mean().expect("Unexpected Empty Row")
            )
        .iter()
        .map(|x| { if *x == 0. { 1. } else { *x }})
        .collect()
}

/// Performs the median ratio normalization method. 
/// Read counts are adjusted by the median of the size factors of each `sgRNA`.
/// Size factors are computed as the geometric mean of the `sgRNAs` across all
/// experimental libraries.
pub fn median_ratio_normalization(
    matrix: &Array2<f64>) -> Array2<f64>
{
    let gmeans = geometric_means(matrix);
    let transformed_matrix = matrix.t().div(gmeans).reversed_axes();
    let sample_scalars = transformed_matrix
        .map_axis(Axis(0), |axis| 1. / median(&axis));
    matrix.mul(sample_scalars)
}


#[cfg(test)]
mod testing {
    use ndarray::{Array1, Array2};
    use ndarray_rand::{RandomExt, rand_distr::Uniform};
    use super::{median, median_ratio_normalization};

    #[test]
    fn test_median_odd() {
        let arr = Array1::range(1., 5., 1.);
        assert_eq!(median(&arr.view()), 3.);
    }

    #[test]
    fn test_median_even() {
        let arr = Array1::range(1., 6., 1.);
        assert_eq!(median(&arr.view()), 3.5);
    }

    #[test]
    fn test_median_normalization() {
        (0..1000).for_each(|_| {
            let matrix = Array2::random((10, 4), Uniform::new(0., 10.));
            let norm = median_ratio_normalization(&matrix);

            // matrices must be equal shape
            assert_eq!(norm.shape(), matrix.shape());
        })
    }
}
