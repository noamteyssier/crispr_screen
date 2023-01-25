use ndarray::{s, Array2, Array1, Axis};
use crate::utils::logging::Logger;
use super::LoggedOls;
use crate::norm::median;

/// Model Mean Variance using Ordinary Least Squares Regression
pub fn model_mean_variance(
    normed_matrix: &Array2<f64>,
    n_controls: usize,
    logger: &Logger) -> Array1<f64>
{
    let model_matrix = if n_controls == 1 {
        normed_matrix.view()
    } else {
        normed_matrix.slice(s![.., ..n_controls])
    };
    let model_mean = model_matrix
        .map_axis(Axis(1), |x| median(&x));
    let model_var = model_matrix.var_axis(Axis(1), 1.);

    let control_mean = normed_matrix
        .slice(s![.., ..n_controls])
        .map_axis(Axis(1), |x| median(&x));
    let logged_ols = LoggedOls::fit(&model_mean, &model_var, logger);
    
    logged_ols.predict(&control_mean)
}

