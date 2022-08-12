use std::fmt::Debug;
use hashbrown::HashSet;
use colored::Colorize;

use crate::norm::Normalization;

pub struct Logger {
    verbose: bool
}
impl Logger {

    pub fn new() -> Self {
        Self { verbose: true }
    }

    fn write_to_stderr<V: Debug>(prompt: &str, value: V) {
        eprintln!(
            "{:width$} {}", 
            format!(">> {}", prompt.bright_green()),
            format!("{:?}", value).white().bold(),
            width = 30
            )
    }

    pub fn num_sgrnas(&self, x: &[String]) {
        if self.verbose {
            Self::write_to_stderr("Number of sgRNAs     : ", x.len());
        }
    }

    pub fn num_genes(&self, x: &[String]) {
        if self.verbose {
            let unique_genes: HashSet<String> = HashSet::from_iter(x.iter().cloned());
            Self::write_to_stderr("Number of Genes      : ", unique_genes.len());
        }
    }

    pub fn norm_method(&self, n: &Normalization) {
        if self.verbose {
            Self::write_to_stderr("Normalization Method : ", n)
        }
    }
}
