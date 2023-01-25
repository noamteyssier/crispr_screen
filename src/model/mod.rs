use clap::ValueEnum;

mod model_mean_variance;
mod logged_ols;
mod ols;
mod wols;

pub use model_mean_variance::model_mean_variance;
use ols::Ols;
use wols::Wols;
use logged_ols::LoggedOls;

#[derive(ValueEnum, Debug, Clone)]
pub enum ModelChoice {
    Ols,
    Wols
}
