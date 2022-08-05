use polars::prelude::{CsvReader, SerReader, DataFrame, PolarsError};
mod math;
mod utils;
mod differential_expression;
use differential_expression::{mageck, Normalization};

fn load_dataframe(path: &str) -> Result<DataFrame, PolarsError>
{
    CsvReader::from_path(path)?
        .has_header(true)
        .with_delimiter(b'\t')
        .finish()
}

fn main() {
    let path = "EXP114.results.tab";
    let frame = load_dataframe(path).unwrap();
    let labels_controls = vec!["EXP114_LOW".to_string()];
    let labels_treatments = vec!["EXP114_HIGH".to_string()];
    mageck(
        &frame,
        &labels_controls,
        &labels_treatments,
        Normalization::MedianRatio
    ).unwrap();
}
