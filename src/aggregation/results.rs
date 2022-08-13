use ndarray::Array1;

pub struct AggregationResult {
    genes: Vec<String>,
    pvalues_low: Array1<f64>,
    pvalues_high: Array1<f64>,
    aggregation_score_low: Array1<f64>,
    aggregation_score_high: Array1<f64>
}
impl AggregationResult {
    pub fn new(
        genes: Vec<String>,
        pvalues_low: Array1<f64>,
        pvalues_high: Array1<f64>,
        aggregation_score_low: Array1<f64>,
        aggregation_score_high: Array1<f64>) -> Self {
        Self {
            genes,
            pvalues_low,
            pvalues_high,
            aggregation_score_low,
            aggregation_score_high,
        }
    }

    pub fn genes(&self) -> &Vec<String> {
        &self.genes
    }

    pub fn pvalues_low(&self) -> &Array1<f64> {
        &self.pvalues_low
    }

    pub fn pvalues_high(&self) -> &Array1<f64> {
        &self.pvalues_high
    }

    pub fn score_low(&self) -> &Array1<f64> {
        &self.aggregation_score_low
    }

    pub fn score_high(&self) -> &Array1<f64> {
        &self.aggregation_score_high
    }
}

