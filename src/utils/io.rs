use polars::prelude::{CsvReader, SerReader, DataFrame, PolarsError};

pub fn load_dataframe(path: &str) -> Result<DataFrame, PolarsError>
{
    CsvReader::from_path(path)?
        .has_header(true)
        .with_delimiter(b'\t')
        .finish()
}
