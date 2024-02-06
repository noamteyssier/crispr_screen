use crate::aggregation::GeneAggregationSelection;
use crate::model::ModelChoice;
use crate::norm::Normalization;
use clap::Parser;

use crate::utils::Adjustment;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Filepath of the input count matrix
    #[arg(short, long)]
    pub input: String,

    /// Labels for Control Samples
    #[arg(short, long, num_args=1.., required=true)]
    pub controls: Vec<String>,

    /// Labels for Treatment Samples
    #[arg(short, long, num_args=1.., required=true)]
    pub treatments: Vec<String>,

    /// Output filename prefix
    ///
    /// sgRNA results will be written to <prefix>.sgrna_results.tsv
    ///
    /// gene results will be written to <prefix>.gene_results.tsv
    ///
    /// hits will be written to <prefix>.hits.tsv
    #[arg(short, long, default_value = "./results")]
    pub output: String,

    /// Count normalization configuration
    ///
    /// If high numbers of zeros are encountered the normalization
    /// method will default to `total` normalization.
    #[arg(short, long, default_value = "median-ratio")]
    pub norm: Normalization,

    /// Gene aggregation configuration
    #[arg(short = 'g', long, default_value = "rra")]
    pub agg: GeneAggregationSelection,

    /// Number of permutations to perform in aRRA
    #[arg(short, long, default_value = "100")]
    pub permutations: usize,

    /// Alpha threshold for aRRA
    #[arg(short, long, default_value = "0.25")]
    pub alpha: f64,

    /// Do not adjust alpha threshold for RRA.
    #[arg(long)]
    pub no_adjust_alpha: bool,

    /// Non-targeting control token
    #[arg(long, default_value = "non-targeting")]
    pub ntc_token: String,

    /// FDR-threshold to use in INC + RRA when thresholding
    #[arg(short = 'F', long, default_value = "0.1")]
    pub fdr: f64,

    /// sgRNA group size of pseudogenes to create for INC
    #[arg(short = 'G', long, default_value = "5")]
    pub inc_group_size: usize,

    /// Calculate FDR threshold using product-score in INC instead of the MWU p-values
    #[arg(long)]
    pub inc_product: bool,

    /// Number of draws to use in INC algorithm
    #[arg(long, default_value = "100")]
    pub n_draws: usize,

    /// Do not write logging information
    #[arg(short, long)]
    pub quiet: bool,

    /// Multiple hypothesis correction method
    #[arg(short = 'f', long, default_value = "bh")]
    pub correction: Adjustment,

    /// Least squares model choice
    #[arg(short, long, default_value = "wols")]
    pub model_choice: ModelChoice,

    /// Set the seed of the run
    #[arg(short, long, default_value = "42")]
    pub seed: u64,

    /// Number of threads to use (defaults to all available)
    #[arg(short = 'T', long)]
    pub threads: Option<usize>,
}
