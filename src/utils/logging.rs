use adjustp::Procedure;
use colored::Colorize;
use hashbrown::HashSet;
use std::fmt::Debug;

use crate::{aggregation::GeneAggregation, io::HitList, model::ModelChoice, norm::Normalization};

#[derive(Default)]
pub struct Logger {
    verbose: bool,
}
impl Logger {
    pub fn new() -> Self {
        Self { verbose: true }
    }

    pub fn new_silent() -> Self {
        Self { verbose: false }
    }

    fn write_to_stderr<V: Debug>(prompt: &str, value: V) {
        eprintln!(
            "{:width$} {}",
            format!(">> {}", prompt.bright_green()),
            format!("{value:.5?}").white().bold(),
            width = 30
        );
    }

    pub fn start_mageck(&self) {
        if self.verbose {
            eprintln!("\n{}", "Run Configuration".bold().underline());
        }
    }

    pub fn num_sgrnas(&self, x: &[String]) {
        if self.verbose {
            Self::write_to_stderr("Number of sgRNAs           : ", x.len());
        }
    }

    pub fn group_names(&self, controls: &[String], treatments: &[String]) {
        if self.verbose {
            Self::write_to_stderr("Control Group              : ", controls);
            Self::write_to_stderr("Treatment Group            : ", treatments);
        }
    }

    pub fn num_genes(&self, x: &[String]) {
        if self.verbose {
            let unique_genes = x.iter().cloned().collect::<HashSet<String>>();
            Self::write_to_stderr("Number of Genes            : ", unique_genes.len());
        }
    }

    pub fn norm_method(&self, n: &Normalization) {
        if self.verbose {
            Self::write_to_stderr("Normalization Method       : ", n);
        }
    }

    pub fn aggregation_method(&self, g: &GeneAggregation) {
        if self.verbose {
            Self::write_to_stderr("Aggregation Method         : ", g);
        }
    }

    pub fn correction(&self, correction: Procedure) {
        if self.verbose {
            Self::write_to_stderr("P-Value Correction Method  : ", correction);
        }
    }

    pub fn start_mean_variance(&self) {
        if self.verbose {
            eprintln!("\n{}", "Modeling Mean Variance".bold().underline());
        }
    }

    pub fn num_outliers(&self, x: usize) {
        if self.verbose {
            Self::write_to_stderr("Removed Outlier sgRNAs     : ", x);
        }
    }

    pub fn num_zeros(&self, x: usize) {
        if self.verbose {
            Self::write_to_stderr("Removed Zero sgRNAs        : ", x);
        }
    }

    pub fn num_varied(&self, x: usize) {
        if self.verbose {
            Self::write_to_stderr("Removed Undervaried sgRNAs : ", x);
        }
    }

    pub fn ols_parameters(&self, model_choice: &ModelChoice, kappa: f64, beta: f64) {
        if self.verbose {
            Self::write_to_stderr("Linear Model Type          : ", model_choice);
            Self::write_to_stderr("Fit Parameter; K           : ", kappa);
            Self::write_to_stderr("Fit Parameter; B           : ", beta);
        }
    }

    pub fn start_gene_aggregation(&self) {
        if self.verbose {
            eprintln!("\n{}", "Performing Gene Aggregation".bold().underline());
        }
    }

    pub fn report_rra_alpha(&self, alpha_low: f64, alpha_high: f64) {
        if self.verbose {
            Self::write_to_stderr("Alpha threshold low        : ", alpha_low);
            Self::write_to_stderr("Alpha threshold high       : ", alpha_high);
        }
    }

    pub fn permutation_sizes(&self, sizes: &[usize]) {
        if self.verbose {
            let mut sorted_sizes = sizes.to_owned();
            sorted_sizes.sort_unstable();
            Self::write_to_stderr("Permutation Sizes          : ", sorted_sizes);
        }
    }

    pub fn report_inc_params(
        &self,
        ntc_token: &str,
        n_genes: usize,
        fdr: f64,
        group_size: usize,
        n_draws: usize,
    ) {
        if self.verbose {
            Self::write_to_stderr("NTC Token                  : ", ntc_token);
            Self::write_to_stderr("Number of Pseudogenes      : ", n_genes);
            Self::write_to_stderr("FDR                        : ", fdr);
            Self::write_to_stderr("Group Size                 : ", group_size);
            Self::write_to_stderr("Number of Draws            : ", n_draws);
        }
    }

    pub fn report_inc_low_threshold(&self, threshold: f64, use_product: bool) {
        if self.verbose {
            if use_product {
                Self::write_to_stderr("Low Product Threshold      : ", threshold);
            } else {
                Self::write_to_stderr("Low Pvalue Threshold       : ", threshold);
            }
        }
    }

    pub fn report_inc_high_threshold(&self, threshold: f64, use_product: bool) {
        if self.verbose {
            if use_product {
                Self::write_to_stderr("High Product Threshold     : ", threshold);
            } else {
                Self::write_to_stderr("High Pvalue Threshold      : ", threshold);
            }
        }
    }

    pub fn report_inc_ntc_std(&self, ntc_std: f64) {
        if self.verbose {
            Self::write_to_stderr("NTC Log2FC Std. Dev.       : ", ntc_std);
        }
    }

    pub fn convert_normalization(&self) {
        if self.verbose {
            eprintln!(
                "\n{}: {}",
                "Warning".bold().yellow(),
                "Numeric instability found in median-ratio normalization. Performing total normalization instead.".bold()
                );
        }
    }

    pub fn hit_list(&self, hit_list: &HitList) {
        if self.verbose {
            eprintln!("\n{}", "Hits".bold().underline());
            Self::write_to_stderr("Number of Hits             : ", hit_list.num_total());
            Self::write_to_stderr("Number Upregulated Hits    : ", hit_list.num_enrichments());
            Self::write_to_stderr("Number Downregulated Hits  : ", hit_list.num_depletions());
        }
    }
}

#[cfg(test)]
mod testing {

    use super::Logger;
    use crate::aggregation::GeneAggregation;
    use crate::model::ModelChoice;
    use crate::norm::Normalization;
    use adjustp::Procedure;

    #[test]
    fn test_logger() {
        let logger = Logger::new();
        logger.start_mageck();
        logger.num_sgrnas(&["a".to_string(), "b".to_string()]);
        logger.num_genes(&["a".to_string(), "b".to_string()]);
        logger.norm_method(&Normalization::MedianRatio);
        logger.aggregation_method(&GeneAggregation::Inc {
            token: "ntc",
            fdr: 0.05,
            group_size: 5,
            use_product: true,
            n_draws: 100,
        });
        logger.correction(Procedure::Bonferroni);
        logger.start_mean_variance();
        logger.num_outliers(1);
        logger.num_zeros(1);
        logger.num_varied(1);
        logger.ols_parameters(&ModelChoice::Ols, 1.0, 1.0);
        logger.start_gene_aggregation();
        logger.report_rra_alpha(1.0, 1.0);
        logger.permutation_sizes(&[1, 2, 3]);
        logger.report_inc_params("NTC", 1, 1.0, 1, 100);
        logger.report_inc_low_threshold(1.0, false);
        logger.report_inc_high_threshold(1.0, false);
        logger.report_inc_low_threshold(1.0, true);
        logger.report_inc_high_threshold(1.0, true);
        logger.convert_normalization();
    }

    #[test]
    fn test_logger_quiet() {
        let logger = Logger::new_silent();
        logger.start_mageck();
        logger.num_sgrnas(&["a".to_string(), "b".to_string()]);
        logger.num_genes(&["a".to_string(), "b".to_string()]);
        logger.norm_method(&Normalization::MedianRatio);
        logger.aggregation_method(&GeneAggregation::Inc {
            token: "ntc",
            fdr: 0.05,
            group_size: 5,
            use_product: false,
            n_draws: 100,
        });
        logger.correction(Procedure::Bonferroni);
        logger.start_mean_variance();
        logger.num_outliers(1);
        logger.num_zeros(1);
        logger.num_varied(1);
        logger.ols_parameters(&ModelChoice::Ols, 1.0, 1.0);
        logger.start_gene_aggregation();
        logger.report_rra_alpha(1.0, 1.0);
        logger.permutation_sizes(&[1, 2, 3]);
        logger.report_inc_params("NTC", 1, 1.0, 1, 100);
        logger.report_inc_low_threshold(1.0, false);
        logger.report_inc_high_threshold(1.0, false);
        logger.report_inc_low_threshold(1.0, true);
        logger.report_inc_high_threshold(1.0, true);
        logger.convert_normalization();
    }
}
