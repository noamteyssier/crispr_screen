use adjustp::{adjust, Procedure};
use ndarray::Array1;
pub struct EnrichmentResult {
    pvalues_low: Array1<f64>,
    pvalues_high: Array1<f64>,
    pvalues_twosided: Array1<f64>,
    fdr: Array1<f64>
}
impl EnrichmentResult {
    pub fn new(pvalues_low: Array1<f64>, pvalues_high: Array1<f64>) -> Self {
        let pvalues_twosided = Self::calculate_twosided(&pvalues_low, &pvalues_high);
        let fdr = Self::calculate_fdr(&pvalues_twosided);
        Self {
            pvalues_low,
            pvalues_high,
            pvalues_twosided,
            fdr
        }
    }

    fn calculate_twosided(pvalues_low: &Array1<f64>, pvalues_high: &Array1<f64>) -> Array1<f64>{
        pvalues_low
            .iter()
            .zip(pvalues_high.iter())
            .map(|(l, h)| l.min(*h) * 2.)
            .collect()
    }

    fn calculate_fdr(pvalues: &Array1<f64>) -> Array1<f64> {
        Array1::from_vec(
            adjust(pvalues.as_slice().unwrap(), Procedure::BenjaminiHochberg)
        )
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
}
