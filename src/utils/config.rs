use adjustp::Procedure;
use crate::{
    norm::Normalization, 
    aggregation::GeneAggregation, 
    model::ModelChoice
};

pub struct Configuration<'a> {
    normalization: Normalization,
    aggregation: GeneAggregation<'a>,
    correction: Procedure,
    model_choice: ModelChoice,
    prefix: &'a str,
}
impl <'a> Configuration<'a> {
    pub fn new(
        normalization: Normalization, 
        aggregation: GeneAggregation<'a>, 
        correction: Procedure, 
        model_choice: ModelChoice,
        prefix: &'a str,
) -> Self {
        Self {
            normalization,
            aggregation,
            correction,
            model_choice,
            prefix,
        }
    }
    pub fn normalization(&self) -> &Normalization {
        &self.normalization
    }
    pub fn aggregation(&self) -> &GeneAggregation {
        &self.aggregation
    }
    pub fn correction(&self) -> Procedure {
        self.correction
    }
    pub fn model_choice(&self) -> &ModelChoice {
        &self.model_choice
    }
    pub fn prefix(&self) -> &str {
        &self.prefix
    }
}
