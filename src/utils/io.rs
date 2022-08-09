use std::fs::File;
use polars::prelude::{CsvReader, CsvWriter, SerReader, SerWriter, DataFrame, PolarsError};

pub fn load_dataframe(path: &str) -> Result<DataFrame, PolarsError>
{
    CsvReader::from_path(path)?
        .has_header(true)
        .with_delimiter(b'\t')
        .finish()
}

fn write_dataframe(name: &str, frame: &mut DataFrame) -> Result<(), PolarsError>
{
    let file = File::create(name)?;
    CsvWriter::new(file)
        .has_header(true)
        .with_delimiter(b'\t')
        .finish(frame)
}

pub fn write_sgrna_results(prefix: &str, frame: &mut DataFrame) -> Result<(), PolarsError>
{
    let filename = format!("{}.sgrna_results.tab", prefix);
    write_dataframe(&filename, frame)
}

pub fn write_gene_results(prefix: &str, frame: &mut DataFrame) -> Result<(), PolarsError>
{
    let filename = format!("{}.gene_results.tab", prefix);
    write_dataframe(&filename, frame)
}

