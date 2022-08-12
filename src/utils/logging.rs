use hashbrown::HashSet;
use colored::Colorize;

pub struct Logger {
    verbose: bool
}
impl Logger {

    pub fn new() -> Self {
        Self { verbose: true }
    }

    fn write_to_stderr(prompt: &str, value: usize) {
        eprintln!(
            "{} {}", 
            format!(">> {}", prompt.bright_green()),
            format!("{}", value).white().bold()
            )
    }

    pub fn num_sgrnas(&self, x: &[String]) {
        if self.verbose {
            Self::write_to_stderr("Number of sgRNAs : ", x.len());
        }
    }

    pub fn num_genes(&self, x: &[String]) {
        if self.verbose {
            let unique_genes: HashSet<String> = HashSet::from_iter(x.iter().cloned());
            Self::write_to_stderr("Number of Genes  : ", unique_genes.len());
        }
    }
}
