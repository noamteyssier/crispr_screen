use anyhow::{bail, Result};
use hashbrown::HashSet;
use ndarray::Array2;
use polars::{
    error::PolarsError,
    frame::DataFrame,
    io::SerReader,
    prelude::{CsvParseOptions, CsvReadOptions},
};
use regex::Regex;
use std::path::PathBuf;

use crate::aggregation::GeneAggregation;

pub fn load_dataframe(path: PathBuf) -> Result<DataFrame, PolarsError> {
    CsvReadOptions::default()
        .with_has_header(true)
        .with_parse_options(CsvParseOptions::default().with_separator(b'\t'))
        .try_into_reader_with_file_path(Some(path))?
        .finish()
}

pub fn match_headers_from_regex_set(
    dataframe: &DataFrame,
    regexes: &[Regex],
) -> Result<Vec<String>> {
    let headers = dataframe.get_column_names();

    let mut set: Vec<String> = regexes
        .iter()
        .flat_map(|re| {
            headers
                .iter()
                .filter(|x| match re.captures(x) {
                    Some(c) => &c.get(0).unwrap().as_str() == *x,
                    None => false,
                })
                .map(|y| y.to_string())
        })
        .fold(HashSet::new(), |mut set, x| {
            set.insert(x);
            set
        })
        .into_iter()
        .collect();
    set.sort_unstable();
    if set.is_empty() {
        let mut s = String::new();
        s.push_str("No matching headers found for the following regex set:\n");
        regexes.iter().for_each(|re| {
            s.push_str(&format!("{re}\n"));
        });
        bail!(s)
    }
    Ok(set)
}

pub fn validate_ntc(sgrna_names: &[String], config: &GeneAggregation) -> Result<()> {
    match config {
        GeneAggregation::Inc {
            token,
            fdr: _,
            group_size: _,
            n_draws: _,
            use_product: _,
        } => {
            if sgrna_names.iter().any(|x| x.contains(token)) {
                Ok(())
            } else {
                bail!("Non-Targeting Token ({token}) not found in any sgrna names - please use RRA or update the provided token.")
            }
        }
        GeneAggregation::GeoPAGG {
            token,
            weight_config: _,
            fdr: _,
            use_product: _,
            zscore_threshold: _,
        } => {
            if let Some(token) = token {
                if sgrna_names.iter().any(|x| x.contains(token)) {
                    Ok(())
                } else {
                    bail!("Non-Targeting Token ({token}) not found in any sgrna names - please use RRA or update the provided token.")
                }
            } else {
                Ok(())
            }
        }
        _ => Ok(()),
    }
}

pub fn get_string_column(dataframe: &DataFrame, idx: usize) -> Vec<String> {
    dataframe
        .select_at_idx(idx)
        .unwrap()
        .iter()
        .map(|x| x.to_string())
        .map(|x| x.replacen('"', "", 2)) // strips quotes
        .collect()
}

pub fn to_ndarray(dataframe: &DataFrame, labels: &[String]) -> Array2<f64> {
    let mut array = Array2::zeros((dataframe.height(), labels.len()));
    for (i, label) in labels.iter().enumerate() {
        let col = dataframe.column(label).unwrap();
        col.f64().unwrap().iter().enumerate().for_each(|(j, x)| {
            array[[j, i]] = x.unwrap_or_default();
        });
    }
    array
}
