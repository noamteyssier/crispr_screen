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
        beta::beta_reg(t_mean.round() + 1.0, r, 1.0 - p)
    } else {
        // implements cumulative density function of the negative binomial distribution without
        // introducing any precision loss (may replace with cdf once updated in statrs)
        beta::beta_reg(r, t_mean.round() + 1.0, p)
    }
}

/// Calculates the negative binomial rate from the provided mean and variance.
/// Calculated using the formula:
/// r = mean**2 / (var - mean)
/// All r that is equal to 0.0 is set to 1.0
fn calculate_r(mean: &Array1<f64>, var: &Array1<f64>) -> Array1<f64> {
    let r = (mean * mean) / (var - mean);
    // r.mapv(|x| if x > 0. {x} else { 1. })
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

/// Performs enrichment testing using a negative binomial distribution
pub fn enrichment_testing(
    normed_matrix: &Array2<f64>,
    adj_var: &Array1<f64>,
    n_controls: usize,
    correction: &Procedure,
) -> EnrichmentResult {
    let control_means = normed_matrix
        .slice(s![.., ..n_controls])
        .map_axis(Axis(1), |x| median(&x));

    let min_control_mean = control_means
        .iter()
        .filter(|x| **x > 0.)
        .copied()
        .reduce(f64::min)
        .expect("Unable to calculate minimum control mean");

    let adj_control_means = control_means.map(|x| if *x == 0. { min_control_mean } else { *x });

    let treatment_means = normed_matrix
        .slice(s![.., n_controls..])
        .map_axis(Axis(1), |x| median(&x));

    let param_r = calculate_r(&adj_control_means, adj_var);
    let param_p = calculate_p(&adj_control_means, adj_var);

    let low = Zip::from(&treatment_means)
        .and(&param_r)
        .and(&param_p)
        .map_collect(|t_mean, r, p| enrichment_test(*t_mean, *r, *p, false));

    let high = Zip::from(&treatment_means)
        .and(&param_r)
        .and(&param_p)
        .map_collect(|t_mean, r, p| enrichment_test(*t_mean, *r, *p, true));

    EnrichmentResult::new(low, high, control_means, treatment_means, correction)
}
