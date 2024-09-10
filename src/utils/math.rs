use anyhow::Result;
use ndarray::{Array1, Array2, Axis};
use rand::Rng;
use rand_distr::{Binomial, Distribution};
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
    let w_sum = weights.sum();
    for (i, row) in m.axis_iter(Axis(0)).enumerate() {
        weighted_sum[i] = ((row.mapv(|x| x.ln()) * weights).sum() / w_sum).exp();
    }
    weighted_sum
}

/// Use `rand_distr` to sample from a binomial distribution
pub fn get_binomial<R: Rng>(probability: f64, n: u64, rng: &mut R) -> Result<u64> {
    let result = Binomial::new(n, probability)?;
    let sample = result.sample(rng);
    Ok(sample)
}

/// Homespun implementation of multinomial distribution using `rand_distr` binomial distribution
/// for use of BTPE algorithm.
///
/// See discussion: https://github.com/statrs-dev/statrs/issues/183#issuecomment-2270083559
pub fn get_multinomial<R: Rng>(pix: &[f64], n: u64, rng: &mut R) -> Result<Vec<u64>> {
    let mut remaining_p = 1.0;
    let mut dn = n;
    let mut mnix = vec![0; pix.len()];
    let size_of_pvec = pix.len();

    for i in 0..size_of_pvec - 1 {
        let binomial = get_binomial(pix[i] / remaining_p, dn, rng)?;
        mnix[i] = binomial;
        dn -= binomial;
        if dn == 0 {
            break;
        }
        remaining_p -= pix[i];
    }
    if dn > 0 {
        mnix[size_of_pvec - 1] = dn;
    }
    Ok(mnix)
}

#[cfg(test)]
mod testing {
    use super::*;
    use ndarray::array;
    use ndarray::Array1;
    use ndarray_rand::{rand_distr::Normal, RandomExt};

    #[test]
    fn test_zscore() {
        let x = Array1::random(100_000, Normal::new(100., 300.).unwrap());
        let z = zscore_transform(&x);

        assert!(z.mean() < Some(1e-6));
        assert!(z.std(0.) - 1. < 1e-6);
    }

    #[test]
    fn test_geometric_mean_weighted() {
        let x = array![[18., 1327., 1024., 1001., 1116.]];
        let w: Array1<f64> = array![0.1, 0.2, 0.3, 0.4, 0.5];
        let gmean = weighted_geometric_mean(&x, &w);
        let expected: Array1<f64> = array![828.1927675627339];
        assert_eq!(gmean, expected);
    }

    #[test]
    fn test_geometric_mean_unweighted() {
        let x = array![[18., 1327., 1024., 1001., 1116.]];
        let w: Array1<f64> = array![1.0, 1.0, 1.0, 1.0, 1.0];
        let gmean = weighted_geometric_mean(&x, &w);
        let expected: Array1<f64> = array![486.7526575476501];
        assert_eq!(gmean, expected);
    }
}
