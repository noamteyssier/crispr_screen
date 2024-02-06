use crate::{
    aggregation::GeneAggregationSelection, model::ModelChoice, norm::Normalization,
    utils::Adjustment,
};
use clap::{Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[clap(subcommand)]
    pub command: Commands,
}

#[derive(Parser, Debug)]
pub struct InputArgs {
    /// Filepath of the input count matrix
    #[arg(short, long)]
    pub input: String,

    /// Labels for Control Samples
    #[arg(short, long, num_args=1.., required=true)]
    pub controls: Vec<String>,

    /// Labels for Treatment Samples
    #[arg(short, long, num_args=1.., required=true)]
    pub treatments: Vec<String>,
}

#[derive(Parser, Debug)]
pub struct DiffAbundanceArgs {
    /// Count normalization configuration
    ///
    /// If high numbers of zeros are encountered the normalization
    /// method will default to `total` normalization.
    #[arg(short, long, default_value = "median-ratio")]
    pub norm: Normalization,

    /// Least squares model choice
    #[arg(short, long, default_value = "wols")]
    pub model_choice: ModelChoice,
}

#[derive(Parser, Debug)]
pub struct RraArgs {
    /// Number of permutations to perform in aRRA
    #[arg(short, long, default_value = "100")]
    pub permutations: usize,

    /// Alpha threshold for aRRA
    #[arg(short, long, default_value = "0.25")]
    pub alpha: f64,

    /// Do not adjust alpha threshold for RRA.
    #[arg(long)]
    pub no_adjust_alpha: bool,

    /// fdr-threshold to use in inc + rra when thresholding
    #[arg(short = 'f', long, default_value = "0.1")]
    pub fdr: f64,
}

#[derive(Parser, Debug)]
pub struct IncArgs {
    /// Non-targeting control token
    #[arg(long, default_value = "non-targeting")]
    pub ntc_token: String,

    /// sgRNA group size of pseudogenes to create for INC
    #[arg(short = 'G', long, default_value = "5")]
    pub inc_group_size: usize,

    /// Calculate FDR threshold using product-score in INC instead of the MWU p-values
    #[arg(long)]
    pub inc_product: bool,

    /// Number of draws to use in INC algorithm
    #[arg(long, default_value = "100")]
    pub n_draws: usize,

    /// fdr-threshold to use in inc + rra when thresholding
    #[arg(short = 'f', long, default_value = "0.1")]
    pub fdr: f64,
}

#[derive(Parser, Debug)]
pub struct MiscArgs {
    /// Do not write logging information
    #[arg(short, long)]
    pub quiet: bool,

    /// Multiple hypothesis correction method
    #[arg(short = 'f', long, default_value = "bh")]
    pub correction: Adjustment,

    /// Set the seed of the run
    #[arg(short, long, default_value = "42")]
    pub seed: u64,

    /// Number of threads to use (defaults to all available)
    #[arg(short = 'T', long)]
    pub threads: Option<usize>,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Perform a differential abundance analysis
    Test {
        #[clap(flatten)]
        input: InputArgs,

        /// Output filename prefix
        ///
        /// sgRNA results will be written to <prefix>.sgrna_results.tsv
        ///
        /// gene results will be written to <prefix>.gene_results.tsv
        ///
        /// hits will be written to <prefix>.hits.tsv
        #[arg(short = 'o', long, default_value = "./results")]
        prefix: String,

        /// Differential abundance arguments
        #[clap(flatten)]
        diff_args: DiffAbundanceArgs,

        /// Gene aggregation configuration
        #[arg(short = 'g', long, default_value = "rra")]
        agg: GeneAggregationSelection,

        /// RRA arguments
        #[clap(flatten)]
        rra: RraArgs,

        /// INC arguments
        #[clap(flatten)]
        inc: IncArgs,

        /// Misc arguments
        #[clap(flatten)]
        misc: MiscArgs,
    },

    Agg {
        #[clap(short, long)]
        input: String,

        /// Output filename prefix
        ///
        /// gene results will be written to <prefix>.gene_results.tsv
        ///
        /// hits will be written to <prefix>.hits.tsv
        #[arg(short = 'o', long, default_value = "./results")]
        prefix: String,

        /// Gene aggregation configuration
        #[arg(short = 'g', long, default_value = "rra")]
        agg: GeneAggregationSelection,

        /// RRA arguments
        #[clap(flatten)]
        rra: RraArgs,

        /// INC arguments
        #[clap(flatten)]
        inc: IncArgs,

        /// Misc arguments
        #[clap(flatten)]
        misc: MiscArgs,
    },
}
