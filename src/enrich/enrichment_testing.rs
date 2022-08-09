use statrs::function::beta;
use ndarray::{s, Array2, Axis, Array1, Zip};

/// Calculates the negative binomial cumulative distribution if measuring depletion otherwise
/// calculates the negative binomial survival function.
fn enrichment_test(
    t_mean: f64,
    r: f64,
    p: f64,
    use_survival: bool) -> f64
{
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
fn calculate_r(
    mean: &Array1<f64>, 
    var: &Array1<f64>) -> Array1<f64>
{
    let r = (mean * mean) / (var - mean);
    r.mapv(|x| if x > 0. {x} else { 1. })
}

/// Calculates the negative binomial probability from the provided mean and variance.
/// Calculated using the formula:
/// p = mean / var
fn calculate_p(
    mean: &Array1<f64>,
    var: &Array1<f64>) -> Array1<f64>
{
    mean / var
}

/// Performs enrichment testing using a negative binomial distribution
pub fn enrichment_testing(
    normed_matrix: &Array2<f64>,
    adj_var: &Array1<f64>,
    n_controls: usize) -> (Array1<f64>, Array1<f64>)
{
    let control_means = normed_matrix
        .slice(s![.., ..n_controls])
        .mean_axis(Axis(1))
        .expect("Unexpected Empty Control Matrix")
        .map(|x| if *x == 0. {1.} else {*x});
    let treatment_means = normed_matrix
        .slice(s![.., n_controls..])
        .mean_axis(Axis(1))
        .expect("Unexpected Empty Treatment Matrix");
    
    let param_r = calculate_r(&control_means, adj_var);
    let param_p = calculate_p(&control_means, adj_var);

    let low = Zip::from(&treatment_means)
        .and(&param_r)
        .and(&param_p)
        .map_collect(|t_mean, r, p| enrichment_test(*t_mean, *r, *p, false));

    let high = Zip::from(&treatment_means)
        .and(&param_r)
        .and(&param_p)
        .map_collect(|t_mean, r, p| enrichment_test(*t_mean, *r, *p, true));

    (low, high)
}
