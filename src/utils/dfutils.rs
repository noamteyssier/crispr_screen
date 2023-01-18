use anyhow::Result;
use polars::prelude::{DataFrame, Float64Type};
use ndarray::Array2;
use thiserror::Error;

#[derive(Error, Debug, PartialEq, Eq)]
pub enum FormatError {
    #[error("Unable to parse sgRNA names expected in first column")]
    MalformedSgrnaNames,
    #[error("Unable to parse gene names expected in second column")]
    MalformedGeneNames,
    #[error("Unable to find sample name(s) in column: {sample:?}")]
    MissingSamples {
        sample: String,
    },
    #[error("Unable to parse samples as numerics")]
    MalformedSamples,
}

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

/// Parse gene names from input dataframe
pub fn parse_genes(
    frame: &DataFrame,
    name: &str) -> Result<Vec<String>>
{
    match parse_to_string_vec(frame, name) {
        Ok(v) => Ok(v),
        Err(_) => return Err(FormatError::MalformedGeneNames.into())
    }
}

/// Parse sgrna names from input dataframe
pub fn parse_sgrna(
    frame: &DataFrame,
    name: &str) -> Result<Vec<String>>
{
    match parse_to_string_vec(frame, name) {
        Ok(v) => Ok(v),
        Err(_) => return Err(FormatError::MalformedSgrnaNames.into())
    }
}

/// Parse subdataframe to a 2D ndarray
pub fn parse_to_ndarray(
    frame: &DataFrame, 
    labels: &[String]) -> Result<Array2<f64>>
{
    match frame.select(labels){
        Ok(v) => match v.to_ndarray::<Float64Type>() {
            Ok(a) => match a.iter().all(|x| x.is_nan()) {
                true => Err(FormatError::MalformedSamples.into()),
                false => Ok(a)
            },
            Err(_) => Err(FormatError::MalformedSamples.into())
        },
        Err(e) => {
            Err(FormatError::MissingSamples{sample: e.to_string()}.into())
        }
    }
}

#[cfg(test)]
mod testing {
    use polars::{df, prelude::*};
    use crate::utils::{parse_sgrna, FormatError};
    use super::{parse_genes, parse_to_ndarray};

    #[test]
    fn test_parse_gene() {
        let frame = df![
            "gene" => ["g1", "g2", "g3"],
            "sgrna" => ["s1", "s2", "s3"]
        ].unwrap();
        let v = parse_genes(&frame, "gene").unwrap();
        assert_eq!(v, vec!["g1", "g2", "g3"]);
    }

    #[test]
    fn test_parse_sgrna() {
        let frame = df![
            "gene" => ["g1", "g2", "g3"],
            "sgrna" => ["s1", "s2", "s3"]
        ].unwrap();
        let v = parse_sgrna(&frame, "sgrna").unwrap();
        assert_eq!(v, vec!["s1", "s2", "s3"]);
    }

    #[test]
    fn test_parse_missing_gene() {
        let frame = df![
            "gene" => ["g1", "g2", "g3"],
            "sgrna" => ["s1", "s2", "s3"]
        ].unwrap();
        let variant = parse_genes(&frame, "not_a_column")
            .unwrap_err()
            .downcast::<FormatError>()
            .unwrap();
        assert_eq!(variant, FormatError::MalformedGeneNames);
    }

    #[test]
    fn test_parse_missing_sgrna() {
        let frame = df![
            "gene" => ["g1", "g2", "g3"],
            "sgrna" => ["s1", "s2", "s3"]
        ].unwrap();
        let variant = parse_sgrna(&frame, "not_a_column")
            .unwrap_err()
            .downcast::<FormatError>()
            .unwrap();
        assert_eq!(variant, FormatError::MalformedSgrnaNames);
    }

    #[test]
    fn test_parse_ndarray() {
        let frame = df![
            "sample_1" => [1, 2, 3],
            "sample_2" => [4, 5, 6],
        ].unwrap();
        parse_to_ndarray(
            &frame, 
            &["sample_1".to_owned(), "sample_2".to_owned()]
            ).unwrap();
    }

    #[test]
    fn test_parse_ndarray_missing_sample() {
        let frame = df![
            "sample_1" => [1, 2, 3],
            "sample_2" => [4, 5, 6],
        ].unwrap(); 
        let variant = parse_to_ndarray(&frame, &["not_a_sample".to_owned()])
            .unwrap_err()
            .downcast::<FormatError>()
            .unwrap();
        assert_eq!(variant, FormatError::MissingSamples { sample: "Not found: not_a_sample".to_string() });
    }

    #[test]
    fn test_parse_ndarray_malformed_array() {
        let frame = df![
            "sample_1" => ["A", "B", "C"],
            "sample_2" => [4, 5, 6],
        ].unwrap(); 
        let variant = parse_to_ndarray(&frame, &["sample_1".to_owned()])
            .unwrap_err()
            .downcast::<FormatError>()
            .unwrap();
        assert_eq!(variant, FormatError::MalformedSamples);
    }
}
