use std::ops::Mul;

use ndarray::Array1;
use ndarray_rand::rand_distr::num_traits::Pow;

/// An implementation of Ordinary Least Squares for a simple linear regression:
/// https://en.wikipedia.org/wiki/Simple_linear_regression
pub struct OLS {
    alpha: f64,
    beta: f64
}
impl OLS {
    /// An implementation of Ordinary Least Squares for a simple linear regression:
    /// https://en.wikipedia.org/wiki/Simple_linear_regression
    pub fn fit(
        x: &Array1<f64>,
        y: &Array1<f64>) -> Self 
    {
        assert_eq!(x.len(), y.len(), "Provided vectors are of unequal size");
        let n = x.len() as f64;

        let mean_x = x.mean().expect("Empty Array in OLS: X");
        let mean_y = y.mean().expect("Empty Array in OLS: Y");
        
        // numerators
        let num1 = x.mul(y).sum();
        let num2 = x.sum();
        let num3 = y.sum();

        // denominators
        let den1 = x.mapv(|x| x.pow(2)).sum();
        let den2 = num2.pow(2);

        // calculate coefficients
        let beta = ((n * num1) - (num2 * num3)) / ((n * den1) - den2);
        let alpha = mean_y - (beta * mean_x);

        Self { alpha, beta }
    }

    /// Predicts the independent variables provided some array of dependent variables
    #[allow(dead_code)]
    pub fn predict(&self, x: &Array1<f64>) -> Array1<f64>
    {
        self.alpha + (self.beta * x)
    }

    /// Calculate the residuals of the model provided the inputs
    #[allow(dead_code)]
    pub fn residuals(&self, x: &Array1<f64>, y: &Array1<f64>) -> Array1<f64>
    {
        y - self.predict(x)
    }

    /// Return Fit Intercept
    pub fn alpha(&self) -> f64
    {
        self.alpha
    }

    /// Return Fit Coefficient
    pub fn beta(&self) -> f64
    {
        self.beta
    }
}


#[cfg(test)]
mod testing {
    use std::ops::Mul;
    use ndarray::Array1;
    use ndarray_rand::{RandomExt, rand_distr::Uniform};
    use super::OLS;
    const EPSILON: f64 = 1e-6;

    #[test]
    pub fn test_ols() {
        let x = Array1::range(0., 100., 1.);
        let y = &x * Array1::random(100, Uniform::new(0., 50.));
        let ols = OLS::fit(&x, &y);
        let res = ols.residuals(&x, &y);

        // the sum of the residuals should be zero
        assert!(res.sum() < EPSILON);

        // the residuals and x values should be uncorrelated
        assert!(x.mul(res).sum() < EPSILON);
    }
}
