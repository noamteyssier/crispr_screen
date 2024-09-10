use crate::{
    aggregation::{GeneAggregationSelection, GeoPAGGWeightConfigEnum},
    enrich::TestStrategy,
    model::ModelChoice,
    norm::Normalization,
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
#[clap(next_help_heading = "Input Arguments")]
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
#[clap(next_help_heading = "Differential Abundance Arguments")]
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

    /// Minimum Base Mean to consider for differential abundance
    #[arg(short = 'M', long, default_value = "100")]
    pub min_base_mean: f64,

    /// Sample testing strategy
    #[arg(short = 'S', long, default_value = "cm")]
    pub strategy: TestStrategy,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "alpha-RRA Arguments")]
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
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "INC Arguments")]
pub struct IncArgs {
    /// sgRNA group size of pseudogenes to create for INC
    #[arg(short = 'G', long, default_value = "5")]
    pub inc_group_size: usize,

    /// Calculate FDR threshold using product-score in INC instead of the MWU p-values
    #[arg(long)]
    pub inc_product: bool,

    /// Number of draws to use in INC algorithm
    #[arg(long, default_value = "100")]
    pub n_draws: usize,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "GeoPAGG Arguments")]
pub struct GeopaggArgs {
    /// Weight configuration for GeoPAGG
    #[arg(long, default_value = "drop-first")]
    pub weight_config: GeoPAGGWeightConfigEnum,

    /// Drop-First weight configuration alpha parameter (only used if weight_config is drop-first)
    #[arg(long, default_value = "0.5")]
    pub df_alpha: f64,

    /// Use all sgRNAs when making amalgam genes
    #[arg(long)]
    pub use_all: bool,

    /// Calculate Empirical FDR using product-score instead of the aggregated p-values
    #[arg(long)]
    pub use_product: bool,

    /// Set a z-score threshold for non-targeting control distribution before making amalgam genes
    #[arg(long)]
    pub zscore_threshold: Option<f64>,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "Miscellaneous Arguments")]
pub struct MiscArgs {
    /// fdr-threshold to use in inc + rra when thresholding
    #[arg(short = 'f', long, default_value = "0.1")]
    pub fdr: f64,

    /// Non-targeting control token
    #[arg(long, default_value = "non-targeting")]
    pub ntc_token: String,

    /// Do not write logging information
    #[arg(short, long)]
    pub quiet: bool,

    /// Multiple hypothesis correction method
    #[arg(short = 'C', long, default_value = "bh")]
    pub correction: Adjustment,

    /// Set the seed of the run
    #[arg(short, long, default_value = "42")]
    pub seed: u64,

    /// Number of threads to use (defaults to all available)
    #[arg(short = 'T', long)]
    pub threads: Option<usize>,
}

#[derive(Parser, Debug)]
#[clap(next_help_heading = "sgRNA Column Name Arguments")]
pub struct SgrnaColumns {
    /// Column name for the low-side p-value
    #[arg(long, default_value = "pvalue_low")]
    pub pvalue_low: String,

    /// Column name for the high-side p-value
    #[arg(long, default_value = "pvalue_high")]
    pub pvalue_high: String,

    /// Column name for the mean value of the base samples
    #[arg(long, default_value = "control")]
    pub control_mean: String,

    /// Column name for the mean value of the treatment samples
    #[arg(long, default_value = "treatment")]
    pub treatment_mean: String,

    /// Column name for the sgRNA names
    #[arg(long, default_value = "sgrna")]
    pub sgrna: String,

    /// Column name for the gene names
    #[arg(long, default_value = "gene")]
    pub gene: String,
}

#[derive(Subcommand, Debug)]
#[clap(next_help_heading = "Subcommands")]
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

        /// GeoPAGG arguments
        #[clap(flatten)]
        geopagg: GeopaggArgs,

        /// Misc arguments
        #[clap(flatten)]
        misc: MiscArgs,

        /// Skip performing gene aggregation
        #[clap(long)]
        skip_agg: bool,
    },

    /// Perform just the gene aggregation given sgRNA results
    Agg {
        /// Filepath of the input sgRNA results
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

        /// GeoPAGG arguments
        #[clap(flatten)]
        geopagg: GeopaggArgs,

        /// Misc arguments
        #[clap(flatten)]
        misc: MiscArgs,

        /// Column names for sgrna results
        #[clap(flatten)]
        columns: SgrnaColumns,
    },

    /// Resample the input count matrix with various parameterizations
    Resample {
        /// Filepath of the input count matrix
        #[clap(short, long)]
        input: String,

        /// Filepath to write the resampled count matrix
        ///
        /// [default: stdout]
        #[arg(short, long)]
        output: Option<String>,

        /// Number of resamples to perform
        #[arg(short, long)]
        n_resamples: usize,

        /// Sequencing Depth to use for resampling
        ///
        /// [default: mean[original_depth]]
        #[arg(short, long)]
        depth: Option<usize>,

        /// Samples to use for calculating mean depth
        ///
        /// Will use the mean depth across the samples provided
        #[arg(short = 'S', long, num_args=1.., conflicts_with = "depth")]
        depth_samples: Option<Vec<String>>,

        /// Seed to use for resampling
        #[arg(long, default_value = "42")]
        seed: u64,

        /// Sample names whose data should be resampled
        #[arg(short, long, num_args=1.., required = true)]
        samples: Vec<String>,
    },
}
