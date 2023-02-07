use adjustp::{adjust, Procedure};
use ndarray::Array1;
pub struct EnrichmentResult {
    pvalues_low: Array1<f64>,
    pvalues_high: Array1<f64>,
    pvalues_twosided: Array1<f64>,
    fdr: Array1<f64>,
    base_means: Array1<f64>,
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
        let base_means = Self::calculate_base_mean(&control_means, &treatment_means);
        let pvalues_twosided = Self::calculate_twosided(&pvalues_low, &pvalues_high);
        let fdr = Self::calculate_fdr(&pvalues_twosided, correction);
        Self {
            pvalues_low,
            pvalues_high,
            pvalues_twosided,
            fdr,
            base_means,
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

    fn calculate_base_mean(control: &Array1<f64>, treatment: &Array1<f64>) -> Array1<f64> {
        (control + treatment) / 2.
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

    pub fn base_means(&self) -> &Array1<f64> {
        &self.base_means
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

#[cfg(test)]
mod testing {
    use adjustp::Procedure;
    use ndarray::{arr1, Zip};
    use ndarray_rand::rand_distr::num_traits::Float;

    #[test]
    fn test_enrichment_result() {
        let pvalues_low = arr1(&[0.1, 0.2, 0.3, 0.4, 0.5]);
        let pvalues_high = arr1(&[0.2, 0.3, 0.4, 0.5, 0.6]);
        let control_means = arr1(&[1., 2., 3., 4., 5.]);
        let treatment_means = arr1(&[2., 3., 4., 5., 6.]);
        let correction = Procedure::BenjaminiHochberg;
        let result = super::EnrichmentResult::new(
            pvalues_low,
            pvalues_high,
            control_means,
            treatment_means,
            correction,
        );
        assert_eq!(result.pvalues_low(), &arr1(&[0.1, 0.2, 0.3, 0.4, 0.5]));
        assert_eq!(result.pvalues_high(), &arr1(&[0.2, 0.3, 0.4, 0.5, 0.6]));
        assert_eq!(result.control_means(), &arr1(&[1., 2., 3., 4., 5.]));
        assert_eq!(result.treatment_means(), &arr1(&[2., 3., 4., 5., 6.]));
    }

    #[test]
    fn test_calculate_twosided() {
        let pvalues_low = arr1(&[0.1, 0.2, 0.3, 0.4, 0.5]);
        let pvalues_high = arr1(&[0.2, 0.3, 0.4, 0.5, 0.6]);
        let twosided = super::EnrichmentResult::calculate_twosided(&pvalues_low, &pvalues_high);
        assert_eq!(twosided, arr1(&[0.2, 0.4, 0.6, 0.8, 1.0]));
    }

    #[test]
    fn test_calculate_fdr() {
        let pvalues = arr1(&[0.1, 0.2, 0.3, 0.4, 0.5]);
        let correction = Procedure::BenjaminiHochberg;
        let fdr = super::EnrichmentResult::calculate_fdr(&pvalues, correction);
        assert_eq!(fdr, arr1(&[0.5, 0.5, 0.5, 0.5, 0.5]));
    }

    #[test]
    fn test_calculate_fold_change() {
        let control = arr1(&[1., 2., 3., 4., 5.]);
        let treatment = arr1(&[2., 3., 4., 5., 6.]);
        let expected = Zip::from(&control)
            .and(&treatment)
            .map_collect(|c, t| (t + 1.) / (c + 1.));
        let fold_change = super::EnrichmentResult::calculate_fold_change(&control, &treatment);
        assert_eq!(fold_change, expected);
    }

    #[test]
    fn test_calculate_log_fold_change() {
        let fold_change = arr1(&[2., 3., 4., 5., 6.]);
        let expected = Zip::from(&fold_change).map_collect(|x| x.log2());
        let log_fold_change = super::EnrichmentResult::calculate_log_fold_change(&fold_change);
        assert_eq!(log_fold_change, expected);
    }

    #[test]
    fn test_calculate_base_mean() {
        let control = arr1(&[1., 2., 3., 4., 5.]);
        let treatment = arr1(&[2., 3., 4., 5., 6.]);
        let expected = Zip::from(&control)
            .and(&treatment)
            .map_collect(|c, t| (c + t) / 2.);
        let base_mean = super::EnrichmentResult::calculate_base_mean(&control, &treatment);
        assert_eq!(base_mean, expected);
    }
}
