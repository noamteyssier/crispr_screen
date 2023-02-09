use clap::ValueEnum;

mod logged_ols;
mod model_mean_variance;
mod ols;
mod sqmean;
mod wols;
mod math;

use logged_ols::LoggedOls;
pub use model_mean_variance::model_mean_variance;
use ols::Ols;
use sqmean::Sqmean;
use wols::Wols;
use math::inverse;

#[derive(ValueEnum, Debug, Clone)]
pub enum ModelChoice {
    /// Ordinary least squares
    Ols,

    /// Weighted ordinary least squares
    Wols,

    /// Squared mean
    Sqmean,
}
