mod enrichment_testing;
mod results;
use clap::ValueEnum;
pub use enrichment_testing::enrichment_testing;
pub use results::EnrichmentResult;

#[derive(Debug, Clone, Copy, ValueEnum, Default)]
pub enum TestStrategy {
    /// Take the median of the counts within a group before testing differential abundance
    #[value(name = "cm")]
    #[default]
    CountMedian,

    /// Test each sample individually and then aggregate the p-values with a geometric mean
    #[value(name = "gm")]
    SampleGeometricMean,

    /// Test each sample individually and then aggregate the p-values with the weighted geometric mean
    /// where the weights are calculated from the magnitude of negative-log p-values of each sample
    #[value(name = "wgm")]
    SampleWeightedGeometricMean,
}
