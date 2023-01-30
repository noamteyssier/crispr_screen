use ndarray::Array1;

pub struct Sqmean {
    kappa: f64,
}
impl Sqmean {
    pub fn fit(log_means: &Array1<f64>, log_variances: &Array1<f64>) -> Self {
        let kappa = (log_variances - (2.0 * log_means)).mean().unwrap();
        Self { kappa }
    }
    pub fn kappa(&self) -> f64 {
        self.kappa
    }
}
