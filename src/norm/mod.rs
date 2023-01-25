mod normalize_counts;
mod median_ratio_norm;
mod total_norm;

pub use normalize_counts::{Normalization, normalize_counts};
use median_ratio_norm::median_ratio_normalization;
use total_norm::total_normalization;
pub use median_ratio_norm::median;
