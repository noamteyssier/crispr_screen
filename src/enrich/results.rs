use adjustp::{adjust, Procedure};
use ndarray::Array1;
pub struct EnrichmentResult {
    pvalues_low: Array1<f64>,
    pvalues_high: Array1<f64>,
    pvalues_twosided: Array1<f64>,
    fdr: Array1<f64>,
    control_means: Array1<f64>,
    treatment_means: Array1<f64>,
    fold_change: Array1<f64>,
    log_fold_change: Array1<f64>,
}
impl EnrichmentResult {
    pub fn new(
        pvalues_low: Array1<f64>,
        pvalues_high: Array1<f64>,
        control_means: Array1<f64>,
        treatment_means: Array1<f64>,
        correction: Procedure,
    ) -> Self {
        let fold_change = Self::calculate_fold_change(&control_means, &treatment_means);
        let log_fold_change = Self::calculate_log_fold_change(&fold_change);
        let pvalues_twosided = Self::calculate_twosided(&pvalues_low, &pvalues_high);
        let fdr = Self::calculate_fdr(&pvalues_twosided, correction);
        Self {
            pvalues_low,
            pvalues_high,
            pvalues_twosided,
            fdr,
            control_means,
            treatment_means,
            fold_change,
            log_fold_change,
        }
    }

    fn calculate_twosided(pvalues_low: &Array1<f64>, pvalues_high: &Array1<f64>) -> Array1<f64> {
        pvalues_low
            .iter()
            .zip(pvalues_high.iter())
            .map(|(l, h)| l.min(*h) * 2.)
            .collect()
    }

    fn calculate_fdr(pvalues: &Array1<f64>, correction: Procedure) -> Array1<f64> {
        Array1::from_vec(adjust(pvalues.as_slice().unwrap(), correction))
    }

    fn calculate_fold_change(control: &Array1<f64>, treatment: &Array1<f64>) -> Array1<f64> {
        (treatment + 1.) / (control + 1.)
    }

    fn calculate_log_fold_change(fold_change: &Array1<f64>) -> Array1<f64> {
        fold_change.map(|x| x.log2())
    }

    pub fn pvalues_low(&self) -> &Array1<f64> {
        &self.pvalues_low
    }

    pub fn pvalues_high(&self) -> &Array1<f64> {
        &self.pvalues_high
    }

    pub fn pvalues_twosided(&self) -> &Array1<f64> {
        &self.pvalues_twosided
    }

    pub fn fdr(&self) -> &Array1<f64> {
        &self.fdr
    }

    pub fn control_means(&self) -> &Array1<f64> {
        &self.control_means
    }

    pub fn treatment_means(&self) -> &Array1<f64> {
        &self.treatment_means
    }

    pub fn fold_change(&self) -> &Array1<f64> {
        &self.fold_change
    }

    pub fn log_fold_change(&self) -> &Array1<f64> {
        &self.log_fold_change
    }
}
