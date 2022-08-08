use polars::prelude::{CsvReader, SerReader, DataFrame, PolarsError};
mod math;
mod utils;
mod differential_expression;
use differential_expression::mageck;
use utils::Normalization;

fn load_dataframe(path: &str) -> Result<DataFrame, PolarsError>
{
    CsvReader::from_path(path)?
        .has_header(true)
        .with_delimiter(b'\t')
        .finish()
}

fn main() {
    let path = "EXP114.results.tab";
    let labels_controls = vec!["EXP114_LOW".to_string()];
    let labels_treatments = vec!["EXP114_HIGH".to_string()];
    let norm_option = "median";

    let norm_method = match norm_option {
        "median" => Normalization::MedianRatio,
        "total" => Normalization::Total,
        _ => panic!("Unexpected normalization method provided: {}", norm_option)
    };

    let frame = load_dataframe(path).unwrap();
    mageck(
        &frame,
        &labels_controls,
        &labels_treatments,
        norm_method
    ).unwrap();
}
