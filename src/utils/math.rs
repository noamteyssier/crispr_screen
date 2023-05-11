use ndarray::Array1;
use std::ops::{Div, Sub};

/// Z-Score Transforms the Provided Array
/// # Arguments
/// * `array` - the array to be transformed
pub fn zscore_transform(array: &Array1<f64>) -> Array1<f64> {
    let mu = array.mean().unwrap_or(0.);
    let sigma = array.std(0.);
    array.sub(mu).div(sigma)
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
