use ndarray::{Array1, Array2, Axis};
use std::ops::{Div, Sub};

/// Z-Score Transforms the Provided Array
/// # Arguments
/// * `array` - the array to be transformed
pub fn zscore_transform(array: &Array1<f64>) -> Array1<f64> {
    let mu = array.mean().unwrap_or(0.);
    let sigma = array.std(0.);
    array.sub(mu).div(sigma)
}

/// Normalizes an array to sum to 1
pub fn normalize(array: &Array1<f64>) -> Array1<f64> {
    let sum = array.sum();
    array.mapv(|x| x / sum)
}

/// Takes the sum of all negative log values along an axis
pub fn negative_log_sum(m: &Array2<f64>, axis: Axis) -> Array1<f64> {
    m.map_axis(axis, |row| row.iter().map(|x| -x.ln()).sum())
}

/// Calculate the weighted geometric mean of a 2D array across the rows of the matrix
pub fn weighted_geometric_mean(m: &Array2<f64>, weights: &Array1<f64>) -> Array1<f64> {
    assert_eq!(m.len_of(Axis(1)), weights.len());
    let mut weighted_sum = Array1::zeros(m.len_of(Axis(0)));
    for (i, row) in m.axis_iter(Axis(0)).enumerate() {
        for (j, val) in row.iter().enumerate() {
            weighted_sum[i] += val.ln() * weights[j];
        }
    }
    weighted_sum.mapv(|x: f64| x.exp())
}

#[cfg(test)]
mod testing {
    use super::zscore_transform;
    use ndarray::Array1;
    use ndarray_rand::{rand_distr::Normal, RandomExt};

    #[test]
    fn test_zscore() {
        let x = Array1::random(100_000, Normal::new(100., 300.).unwrap());
        let z = zscore_transform(&x);

        assert!(z.mean() < Some(1e-6));
        assert!(z.std(0.) - 1. < 1e-6);
    }
}
