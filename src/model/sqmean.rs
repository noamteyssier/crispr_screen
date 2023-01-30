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

#[cfg(test)]
mod testing {
    use super::*;
    use ndarray::array;

    #[test]
    fn test_sqmean() {
        let log_means = array![1.0, 2.0, 3.0];
        let log_variances = array![3.0, 5.0, 7.0];
        let sqmean = Sqmean::fit(&log_means, &log_variances);
        assert_eq!(sqmean.kappa(), 1.0);
    }

}
