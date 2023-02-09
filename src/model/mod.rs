use clap::ValueEnum;

mod logged_ols;
mod math;
mod model_mean_variance;
mod ols;
mod sqmean;
mod wols;

use logged_ols::LoggedOls;
use math::inverse;
pub use model_mean_variance::model_mean_variance;
use ols::Ols;
use sqmean::Sqmean;
use wols::Wols;

#[derive(ValueEnum, Debug, Clone)]
pub enum ModelChoice {
    /// Ordinary least squares
    Ols,

    /// Weighted ordinary least squares
    Wols,

    /// Squared mean
    Sqmean,
}
