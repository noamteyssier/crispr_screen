use adjustp::Procedure;
use colored::Colorize;
use hashbrown::HashSet;
use std::fmt::Debug;

use crate::{aggregation::GeneAggregation, model::ModelChoice, norm::Normalization};

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

    pub fn report_inc_params(&self, ntc_token: &str, n_genes: usize, fdr: f64, group_size: usize) {
        if self.verbose {
            Self::write_to_stderr("NTC Token                  : ", ntc_token);
            Self::write_to_stderr("Number of Pseudogenes      : ", n_genes);
            Self::write_to_stderr("FDR                        : ", fdr);
            Self::write_to_stderr("Group Size                 : ", group_size);
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
}
