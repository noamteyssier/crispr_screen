use std::fmt::Debug;
use hashbrown::HashSet;
use colored::Colorize;

use crate::{norm::Normalization, aggregation::GeneAggregation};

pub struct Logger {
    verbose: bool
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
            format!("{:.3?}", value).white().bold(),
            width = 30
            )
    }

    pub fn start_mageck(&self) {
        if self.verbose {
            eprintln!(
                "\n{}",
                "Run Configuration".bold().underline()
                )
        }
    }

    pub fn num_sgrnas(&self, x: &[String]) {
        if self.verbose {
            Self::write_to_stderr(
                "Number of sgRNAs           : ", 
                x.len());
        }
    }

    pub fn num_genes(&self, x: &[String]) {
        if self.verbose {
            let unique_genes: HashSet<String> = HashSet::from_iter(x.iter().cloned());
            Self::write_to_stderr(
                "Number of Genes            : ", 
                unique_genes.len());
        }
    }

    pub fn norm_method(&self, n: &Normalization) {
        if self.verbose {
            Self::write_to_stderr(
                "Normalization Method       : ", 
                n)
        }
    }

    pub fn aggregation_method(&self, g: &GeneAggregation) {
        if self.verbose {
            Self::write_to_stderr(
                "Aggregation Method         : ", 
                g)
        }
    }

    pub fn start_mean_variance(&self) {
        if self.verbose {
            eprintln!(
                "\n{}",
                "Modeling Mean Variance".bold().underline()
                )
        }
    }

    pub fn num_outliers(&self, x: usize) {
        if self.verbose {
            Self::write_to_stderr(
                "Removed Outlier sgRNAs     : ", 
                x);
        }
    }

    pub fn num_zeros(&self, x: usize) {
        if self.verbose {
            Self::write_to_stderr(
                "Removed Zero sgRNAs        : ", 
                x);
        }
    }

    pub fn num_varied(&self, x: usize) {
        if self.verbose {
            Self::write_to_stderr(
                "Removed Undervaried sgRNAs : ", 
                x);
        }
    }

    pub fn ols_parameters(&self, kappa: f64, beta: f64) {
        if self.verbose {
            Self::write_to_stderr(
                "Fit Parameter; K           : ", 
                kappa);
            Self::write_to_stderr(
                "Fit Parameter; B           : ", 
                beta);
        }
    }

    pub fn start_gene_aggregation(&self) {
        if self.verbose {
            eprintln!(
                "\n{}",
                "Performing Gene Aggregation".bold().underline()
                )
        }
    }

    pub fn permutation_sizes(&self, sizes: &Vec<usize>) {
        if self.verbose {
            let mut sorted_sizes = sizes.clone();
            sorted_sizes.sort_unstable();
            Self::write_to_stderr(
                "Permutation Sizes          : ", 
                sorted_sizes);
        }
    }

    pub fn num_ntcs(&self, num_ntc: usize) {
        if self.verbose {
            Self::write_to_stderr(
                "Number of Found Controls   : ",
                num_ntc
                );
        }
    }

}
