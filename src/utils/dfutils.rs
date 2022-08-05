use anyhow::Result;
use polars::prelude::{DataFrame, Float64Type};
use ndarray::Array2;

/// Parse a series expected to be a container of strings to a vector of strings
pub fn parse_to_string_vec(
    frame: &DataFrame, 
    name: &str) -> Result<Vec<String>>
{
    let vec = frame.column(name)?
        .utf8()?
        .into_iter()
        .map(|x| x.unwrap().to_owned())
        .collect();
    Ok(vec)
}

/// Parse subdataframe to a 2D ndarray
pub fn parse_to_ndarray(
    frame: &DataFrame, 
    labels: &[String]) -> Result<Array2<f64>>
{
    Ok(frame.select(labels)?.to_ndarray::<Float64Type>()?)
}
