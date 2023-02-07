use super::EnrichmentResult;
use crate::norm::median;
use adjustp::Procedure;
use ndarray::{s, Array1, Array2, Axis, Zip};
use statrs::function::beta;

/// Calculates the negative binomial cumulative distribution if measuring depletion otherwise
/// calculates the negative binomial survival function.
fn enrichment_test(t_mean: f64, r: f64, p: f64, use_survival: bool) -> f64 {
    if use_survival {
        // implements survival function of the negative binomial distribution
        beta::beta_reg(t_mean + 1.0, r, 1.0 - p)
    } else {
        // implements cumulative density function of the negative binomial distribution without
        // introducing any precision loss (may replace with cdf once updated in statrs)
        beta::beta_reg(r, t_mean + 1.0, p)
    }
}

/// Calculates the negative binomial rate from the provided mean and variance.
/// Calculated using the formula:
/// r = mean**2 / (var - mean)
/// All r that is equal to 0.0 is set to 1.0
fn calculate_r(mean: &Array1<f64>, var: &Array1<f64>) -> Array1<f64> {
    let r = (mean * mean) / (var - mean);
    r.mapv(|x| if x >= 1. { x } else { 1. })
}

/// Calculates the negative binomial probability from the provided mean and variance.
/// Calculated using the formula:
/// p = mean / var
fn calculate_p(mean: &Array1<f64>, var: &Array1<f64>) -> Array1<f64> {
    (mean / var)
        .into_iter()
        .map(|x| if x <= 1.0 { x } else { 1.0 })
        .collect()
}

/// Calculates the minimum value in an array that is greater than 0.0
fn get_nonzero_minimum(mean: &Array1<f64>) -> f64 {
    mean.iter()
        .filter(|x| **x > 0.)
        .copied()
        .reduce(f64::min)
        .expect("Unable to calculate minimum control mean")
}

/// Sets all values in an array equal to 0.0 to the provided minimum value
fn set_zero_to_minimum(array: &Array1<f64>, minimum: f64) -> Array1<f64> {
    array.map(|x| if *x == 0. { minimum } else { *x })
}

/// Calculates the median of each row in an array
fn row_median(array: &Array2<f64>) -> Array1<f64> {
    array.map_axis(Axis(1), |x| median(&x))
}

/// Selects the first `n_controls` columns from an array
fn select_controls(array: &Array2<f64>, n_controls: usize) -> Array2<f64> {
    array.slice(s![.., ..n_controls]).to_owned()
}

/// Selects the last `n_treatments` columns from an array
fn select_treatments(array: &Array2<f64>, n_controls: usize) -> Array2<f64> {
    array.slice(s![.., n_controls..]).to_owned()
}

/// Maps the enrichment test function over the provided arrays
fn map_enrichment(
    treatment_means: &Array1<f64>,
    param_r: &Array1<f64>,
    param_p: &Array1<f64>,
    survival: bool,
) -> Array1<f64> {
    Zip::from(treatment_means)
        .and(param_r)
        .and(param_p)
        .map_collect(|t_mean, r, p| enrichment_test(*t_mean, *r, *p, survival))
}

/// Performs enrichment testing using a negative binomial distribution
pub fn enrichment_testing(
    normed_matrix: &Array2<f64>,
    adj_var: &Array1<f64>,
    n_controls: usize,
    correction: Procedure,
) -> EnrichmentResult {
    // Subset the control and treatments and calculate the median of each sgrna
    let control_means = row_median(&select_controls(normed_matrix, n_controls));
    let treatment_means = row_median(&select_treatments(normed_matrix, n_controls));

    // Adjust the control means to ensure that all values are greater than 0.0
    let min_control_mean = get_nonzero_minimum(&control_means);
    let adj_control_means = set_zero_to_minimum(&control_means, min_control_mean);

    // Calculate the negative binomial parameters
    let param_r = calculate_r(&adj_control_means, adj_var);
    let param_p = calculate_p(&adj_control_means, adj_var);

    // Perform the enrichment test
    let low = map_enrichment(&treatment_means, &param_r, &param_p, false);
    let high = map_enrichment(&treatment_means, &param_r, &param_p, true);

    EnrichmentResult::new(low, high, control_means, treatment_means, correction)
}

#[cfg(test)]
mod testing {
    use super::enrichment_test;

    fn test_almost<F>(t: f64, r: f64, p: f64, survival: bool, expected: f64, acc: f64, eval: F)
    where
        F: Fn(f64, f64, f64, bool) -> f64,
    {
        let result = eval(t, r, p, survival);
        assert!(
            (result - expected).abs() < acc,
            "Expected {}, got {}",
            expected,
            result
        )
    }

    #[test]
    /// adapted from [statrs negative binomial cdf testing](https://github.com/statrs-dev/statrs/blob/master/src/distribution/negative_binomial.rs#L461)
    fn test_enrichment_test_cdf() {
        test_almost(0.0, 1.0, 0.3, false, 0.3, 1e-08, enrichment_test);
        test_almost(1.0, 1.0, 0.3, false, 0.51, 1e-08, enrichment_test);
        test_almost(4.0, 1.0, 0.3, false, 0.83193, 1e-08, enrichment_test);
    }

    #[test]
    /// adapted from [statrs negative binomial sf testing](https://github.com/statrs-dev/statrs/blob/master/src/distribution/negative_binomial.rs#L474)
    fn test_enrichment_test_sf() {
        test_almost(0.0, 1.0, 0.3, true, 0.7, 1e-08, enrichment_test);
        test_almost(1.0, 1.0, 0.3, true, 0.49, 1e-08, enrichment_test);
        test_almost(
            4.0,
            1.0,
            0.3,
            true,
            0.1680699999999986,
            1e-08,
            enrichment_test,
        );
    }

    #[test]
    fn test_calculate_p() {
        let mean = ndarray::arr1(&[1., 2., 3.]);
        let var = ndarray::arr1(&[2., 4., 6.]);
        let expected = ndarray::arr1(&[0.5, 0.5, 0.5]);
        let result = super::calculate_p(&mean, &var);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_calculate_r() {
        let mean = ndarray::arr1(&[1., 2., 3.]);
        let var = ndarray::arr1(&[2., 4., 6.]);
        let expected = ndarray::arr1(&[1., 2., 3.]);
        let result = super::calculate_r(&mean, &var);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_get_nonzero_minimum() {
        let mean = ndarray::arr1(&[0., 1., 2., 3.]);
        let expected = 1.;
        let result = super::get_nonzero_minimum(&mean);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_set_nonzero_minimum() {
        let mean = ndarray::arr1(&[0., 1., 2., 3.]);
        let expected = ndarray::arr1(&[1., 1., 2., 3.]);
        let result = super::set_zero_to_minimum(&mean, 1.);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_row_median() {
        let mean = ndarray::arr2(&[[0., 1., 2., 3.], [0., 1., 2., 3.]]);
        let expected = ndarray::arr1(&[1.5, 1.5]);
        let result = super::row_median(&mean);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_select_controls() {
        let mean = ndarray::arr2(&[[0., 1., 2., 3.], [0., 1., 2., 3.]]);
        let expected = ndarray::arr2(&[[0., 1.], [0., 1.]]);
        let result = super::select_controls(&mean, 2);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_select_treatments() {
        let mean = ndarray::arr2(&[[0., 1., 2., 3.], [0., 1., 2., 3.]]);
        let expected = ndarray::arr2(&[[2., 3.], [2., 3.]]);
        let result = super::select_treatments(&mean, 2);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_map_enrichment() {
        let mean = ndarray::arr1(&[0., 1., 4.]);
        let r = ndarray::arr1(&[1., 1., 1.]);
        let p = ndarray::arr1(&[0.3, 0.3, 0.3]);
        let expected = ndarray::arr1(&[0.3, 0.51, 0.83193]);
        let result = super::map_enrichment(&mean, &r, &p, false);
        expected
            .iter()
            .zip(result.iter())
            .for_each(|(e, r)| assert!((e - r).abs() < 1e-08));
    }
}
