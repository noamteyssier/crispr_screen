use ndarray::{Array1, Array2, Axis};
use super::inverse;

/// An implementation of [Ordinary Least Squares](https://en.wikipedia.org/wiki/Ordinary_least_squares#Matrix/vector_formulation) using a Matrix/Vector Formulation
#[derive(Debug)]
pub struct Ols {
    alpha: f64,
    beta: f64,
}
impl Ols {
    /// Fits an Ordinary Least Squares Linear Regression
    pub fn fit(x: &Array1<f64>, y: &Array1<f64>) -> Self {
        assert_eq!(
            x.len(),
            y.len(),
            "Provided vectors to OLS are of unequal size"
        );
        let n = x.len();

        // create design matrix with intercept as first column
        let mut mat_x = Array2::ones((n, 1));
        mat_x
            .push_column(x.view())
            .expect("Unable to construct design matrix in OLS");

        // convert y to 2D matrix
        let mat_y = y.clone().insert_axis(Axis(1));

        // B = inv(XtX)XtY
        let xt = mat_x.t();
        let xtx = xt.dot(&mat_x);
        let inv_xtx = inverse(&xtx);

        // drop an axis from the solution for a 1D vector
        let solution = inv_xtx.dot(&xt).dot(&mat_y).remove_axis(Axis(1));

        let alpha = solution[0];
        let beta = solution[1];

        Self { alpha, beta }
    }

    /// Predicts the independent variables provided some array of dependent variables
    #[allow(dead_code)]
    pub fn predict(&self, x: &Array1<f64>) -> Array1<f64> {
        self.alpha + (self.beta * x)
    }

    /// Calculate the residuals of the model provided the inputs
    #[allow(dead_code)]
    pub fn residuals(&self, x: &Array1<f64>, y: &Array1<f64>) -> Array1<f64> {
        y - self.predict(x)
    }

    /// Return Fit Intercept
    pub fn alpha(&self) -> f64 {
        self.alpha
    }

    /// Return Fit Coefficient
    pub fn beta(&self) -> f64 {
        self.beta
    }
}

#[cfg(test)]
mod testing {
    use super::Ols;
    use ndarray::Array1;
    use ndarray_rand::{rand_distr::Uniform, RandomExt};
    use std::ops::{Add, Mul, Sub};
    const EPSILON: f64 = 1e-6;

    #[test]
    pub fn test_ols() {
        let x = Array1::range(0., 100., 1.);
        let y = &x * Array1::random(100, Uniform::new(0., 50.));
        let ols = Ols::fit(&x, &y);
        let res = ols.residuals(&x, &y);

        // the sum of the residuals should be zero
        assert!(res.sum() < EPSILON);

        // the residuals and x values should be uncorrelated
        assert!(x.mul(res).sum() < EPSILON);
    }

    fn f(x: &Array1<f64>, m: f64, b: f64) -> Array1<f64> {
        x.mul(m).add(b)
    }

    #[test]
    pub fn test_ols_known_function() {
        let n = 5000;
        let x = Array1::linspace(0., 10., n);

        for m in Array1::random(5, Uniform::new(5.0, 10.0)) {
            for b in Array1::random(5, Uniform::new(5.0, 10.0)) {
                let y = f(&x, m, b);
                let ols = Ols::fit(&x, &y);
                assert!(ols.alpha().sub(b).abs() < EPSILON);
                assert!(ols.beta().sub(m).abs() < EPSILON);
            }
        }
    }
}
