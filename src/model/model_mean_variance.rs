use ndarray::{s, Array2, Array1, Axis};
use super::LoggedOLS;

/// Model Mean Variance using Ordinary Least Squares Regression
pub fn model_mean_variance(
    normed_matrix: &Array2<f64>,
    n_controls: usize) -> Array1<f64>
{
    let model_matrix = if n_controls == 1 {
        normed_matrix.view()
    } else {
        normed_matrix.slice(s![.., ..n_controls])
    };
    let model_mean = model_matrix.mean_axis(Axis(1))
        .expect("Unexpected Empty Model Matrix");
    let model_var = model_matrix.var_axis(Axis(1), 1.);

    let control_mean = normed_matrix
        .slice(s![.., ..n_controls])
        .mean_axis(Axis(1))
        .expect("Unexpected empty control matrix");
    let logged_ols = LoggedOLS::fit(&model_mean, &model_var);
    
    logged_ols.predict(&control_mean)
}

