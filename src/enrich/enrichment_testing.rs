use super::{EnrichmentResult, TestStrategy};
use crate::{
    norm::median,
    utils::logging::Logger,
    utils::math::{negative_log_sum, normalize, weighted_geometric_mean},
};
use adjustp::Procedure;
use ndarray::{s, stack, Array1, Array2, Axis, Zip};
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

/// Sets all values in an array equal to 0.0 to the minimum nonzero value in the array
fn set_zero_to_minimum_nonzero(array: &Array1<f64>) -> Array1<f64> {
    let minimum = get_nonzero_minimum(array);
    set_zero_to_minimum(array, minimum)
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
    treatment_arr: &Array1<f64>,
    param_r: &Array1<f64>,
    param_p: &Array1<f64>,
    survival: bool,
) -> Array1<f64> {
    Zip::from(treatment_arr)
        .and(param_r)
        .and(param_p)
        .map_collect(|val, r, p| enrichment_test(*val, *r, *p, survival))
}

fn map_enrichment_2d(
    normed_matrix: &Array2<f64>,
    param_p: &Array1<f64>,
    param_r: &Array1<f64>,
    survival: bool,
    weighted: bool,
    logger: &Logger,
) -> Array1<f64> {
    // Map the enrichment test function over each sample independently
    let arrays = normed_matrix
        .axis_iter(Axis(1))
        .map(|col| map_enrichment(&col.to_owned(), param_r, param_p, survival))
        .collect::<Vec<_>>();

    // Stack the results
    let array_views = arrays.iter().map(|x| x.view()).collect::<Vec<_>>();
    let stack = stack(Axis(1), &array_views).unwrap();

    // Calculate the weight of each sample as the sum of the negative log p-values
    // across all treatments
    let weights = if weighted {
        normalize(&negative_log_sum(&stack, Axis(0)))
    } else {
        Array1::ones(stack.len_of(Axis(1)))
    };
    logger.sample_weights(survival, &weights);

    // Perform a weighted geometric mean of the p-values
    // weighted by the magnitude of the negative log p-values
    // across all treatments
    weighted_geometric_mean(&stack, &weights)
}

/// Performs enrichment testing on each sample independently then aggregates the results.
/// Aggregation is performed by first calculating the sum of log p-values for each sample
/// as each samples weight, then performing a weighted geometric mean of the p-values.
pub fn geometric_enrichment_testing(
    normed_matrix: &Array2<f64>,
    adj_var: &Array1<f64>,
    n_controls: usize,
    correction: Procedure,
    weighted: bool,
    logger: &Logger,
) -> EnrichmentResult {
    let treatment_2d = select_treatments(normed_matrix, n_controls);
    let control_means = row_median(&select_controls(normed_matrix, n_controls));
    let treatment_means = row_median(&treatment_2d);
    let min_control_mean = get_nonzero_minimum(&control_means);
    let adj_control_means = set_zero_to_minimum(&control_means, min_control_mean);

    let param_r = calculate_r(&adj_control_means, adj_var);
    let param_p = calculate_p(&adj_control_means, adj_var);

    let low_geom_mean =
        map_enrichment_2d(&treatment_2d, &param_p, &param_r, false, weighted, logger);
    let high_geom_mean =
        map_enrichment_2d(&treatment_2d, &param_p, &param_r, true, weighted, logger);

    EnrichmentResult::new(
        low_geom_mean,
        high_geom_mean,
        control_means,
        treatment_means,
        correction,
    )
}

pub fn median_enrichment_testing(
    normed_matrix: &Array2<f64>,
    adj_var: &Array1<f64>,
    n_controls: usize,
    correction: Procedure,
) -> EnrichmentResult {
    // Subset the control and treatments and calculate the median of each sgrna
    let control_means = row_median(&select_controls(normed_matrix, n_controls));
    let treatment_means = row_median(&select_treatments(normed_matrix, n_controls));

    // Adjust the control means to ensure that all values are greater than 0.0
    let adj_control_means = set_zero_to_minimum_nonzero(&control_means);

    // Calculate the negative binomial parameters
    let param_r = calculate_r(&adj_control_means, adj_var);
    let param_p = calculate_p(&adj_control_means, adj_var);

    // Perform the enrichment test
    let low = map_enrichment(&treatment_means, &param_r, &param_p, false);
    let high = map_enrichment(&treatment_means, &param_r, &param_p, true);

    // Adjust p-values to set zeros to the minimum non-zero value
    let low = set_zero_to_minimum_nonzero(&low);
    let high = set_zero_to_minimum_nonzero(&high);

    EnrichmentResult::new(low, high, control_means, treatment_means, correction)
}

/// Performs enrichment testing using a negative binomial distribution
///
/// Samples are first split into control and treatment groups, then the median of each sgRNA
/// is calculated for each group.
pub fn enrichment_testing(
    normed_matrix: &Array2<f64>,
    adj_var: &Array1<f64>,
    n_controls: usize,
    correction: Procedure,
    strategy: TestStrategy,
    logger: &Logger,
) -> EnrichmentResult {
    logger.start_differential_abundance();
    logger.sample_aggregation_strategy(strategy);
    match strategy {
        TestStrategy::SampleWeightedGeometricMean => geometric_enrichment_testing(
            normed_matrix,
            adj_var,
            n_controls,
            correction,
            true,
            logger,
        ),
        TestStrategy::SampleGeometricMean => geometric_enrichment_testing(
            normed_matrix,
            adj_var,
            n_controls,
            correction,
            false,
            logger,
        ),
        TestStrategy::CountMedian => {
            median_enrichment_testing(normed_matrix, adj_var, n_controls, correction)
        }
    }
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

    #[test]
    fn test_set_zero_to_minimum_nonzero() {
        let x = ndarray::arr1(&[1., 2., 3., 0., 0., 1.]);
        let expected = ndarray::arr1(&[1., 2., 3., 1., 1., 1.]);
        let result = super::set_zero_to_minimum_nonzero(&x);
        assert_eq!(result, expected);
    }
}
