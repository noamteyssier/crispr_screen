use clap::ValueEnum;

mod logged_ols;
mod model_mean_variance;
mod ols;
mod wols;

use logged_ols::LoggedOls;
pub use model_mean_variance::model_mean_variance;
use ols::Ols;
use wols::Wols;

#[derive(ValueEnum, Debug, Clone)]
pub enum ModelChoice {

    /// Ordinary least squares
    Ols,
    
    /// Weighted ordinary least squares
    Wols,
}
