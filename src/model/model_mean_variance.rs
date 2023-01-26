use super::{LoggedOls, ModelChoice};
use crate::norm::median;
use crate::utils::logging::Logger;
use ndarray::{s, Array1, Array2, Axis};

/// Model Mean Variance using Ordinary Least Squares Regression
pub fn model_mean_variance(
    normed_matrix: &Array2<f64>,
    n_controls: usize,
    model_choice: &ModelChoice,
    logger: &Logger,
) -> Array1<f64> {
    let model_matrix = if n_controls == 1 {
        normed_matrix.view()
    } else {
        normed_matrix.slice(s![.., ..n_controls])
    };
    let model_mean = model_matrix.map_axis(Axis(1), |x| median(&x));
    let model_var = model_matrix.var_axis(Axis(1), 1.);

    let control_mean = normed_matrix
        .slice(s![.., ..n_controls])
        .map_axis(Axis(1), |x| median(&x));
    let logged_ols = LoggedOls::fit(&model_mean, &model_var, model_choice, logger);

    logged_ols.predict(&control_mean)
}
