use statrs::function::beta;
use ndarray::{s, Array2, Axis, Array1, Zip};
use ndarray_rand::rand_distr::num_traits::Pow;
use std::ops::{Div, Sub};

/// Calculates the negative binomial cumulative distribution if measuring depletion otherwise
/// calculates the negative binomial survival function.
fn enrichment_test(
    t_mean: f64,
    r: f64,
    p: f64,
    is_positive: bool) -> f64
{
    if is_positive {
        // implements survival function of the negative binomial distribution
        beta::beta_reg(t_mean.round() + 1.0, r, 1.0 - p)
    } else {
        // implements cumulative density function of the negative binomial distribution without
        // introducing any precision loss (may replace with cdf once updated in statrs)
        beta::beta_reg(r, t_mean.round() + 1.0, p)
    }
}

/// Performs enrichment testing using a negative binomial distribution
pub fn enrichment_testing(
    normed_matrix: &Array2<f64>,
    adj_var: &Array1<f64>,
    n_controls: usize) -> Array1<f64>
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
    
    let param_r = (&control_means.mapv(|x| x.pow(2))).div(adj_var.sub(&control_means)).mapv(|x| if x > 0. { x } else { 1. });
    let param_p = (&control_means).div(adj_var);
    let positive_lfc: Array1<bool> = (&control_means).iter().zip(treatment_means.iter()).map(|(c, t)| c < t).collect();

    Zip::from(&treatment_means)
        .and(&param_r)
        .and(&param_p)
        .and(&positive_lfc)
        .map_collect(|t_mean, r, p, is_positive| enrichment_test(*t_mean, *r, *p, *is_positive))
}

