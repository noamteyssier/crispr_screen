use ndarray::{Array1, Array2, Axis};
use ndarray_linalg::solve::Inverse;

/// An implementation of [Weighted Least Squares](https://en.wikipedia.org/wiki/Weighted_least_squares) using a Matrix/Vector Formulation
#[derive(Debug)]
pub struct Wols {
    alpha: f64,
    beta: f64,
}
impl Wols {
    /// Fits a Weighted Least Squares Linear Regression
    pub fn fit(x: &Array1<f64>, y: &Array1<f64>, w: &Array1<f64>) -> Self {
        assert_eq!(
            x.len(),
            y.len(),
            "Provided vectors to WLS are of unequal size"
        );
        assert_eq!(
            x.len(),
            w.len(),
            "Provided weights to WLS are of unequal size"
        );
        let n = x.len();

        // create design matrix with intercept as first column
        let mut mat_x = Array2::ones((n, 1));
        mat_x
            .push_column(x.view())
            .expect("Unable to construct design matrix in OLS");

        // convert y to 2D matrix
        let mat_y = y.clone().insert_axis(Axis(1));

        // create diagonal weight matrix
        let mat_w = Array2::from_diag(w);

        // B = inv(XtWX)XtWY
        let xt = mat_x.t();
        let xtw = xt.dot(&mat_w);
        let xtwx = xtw.dot(&mat_x);
        let inv_xtwx = Inverse::inv(&xtwx).expect("Unable to computer inverse matrix in OLS");

        // drop an axis from the solution for a 1D vector
        let solution = inv_xtwx.dot(&xtw).dot(&mat_y).remove_axis(Axis(1));

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
    use super::Wols;
    use ndarray::Array1;
    use ndarray_rand::{rand_distr::Uniform, RandomExt};
    use std::ops::{Add, Mul, Sub};
    const EPSILON: f64 = 1e-6;

    fn f(x: &Array1<f64>, m: f64, b: f64) -> Array1<f64> {
        x.mul(m).add(b)
    }

    #[test]
    pub fn test_wls_known_function() {
        let n = 5000;
        let x = Array1::linspace(0., 10., n);
        let w = x.clone();

        for m in Array1::random(5, Uniform::new(5.0, 10.0)) {
            for b in Array1::random(5, Uniform::new(5.0, 10.0)) {
                let y = f(&x, m, b);
                let ols = Wols::fit(&x, &y, &w);
                assert!(ols.alpha().sub(b).abs() < EPSILON);
                assert!(ols.beta().sub(m).abs() < EPSILON);
            }
        }
    }
}
