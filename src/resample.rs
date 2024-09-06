use std::{fs::File, io::BufWriter};

use anyhow::Result;
use bon::builder;
use ndarray::prelude::*;
use ndarray_rand::rand::SeedableRng;
use polars::prelude::*;
use rand_chacha::ChaCha8Rng;

use crate::{
    io::load_dataframe,
    norm::{normalize_counts, Normalization},
    utils::{logging::Logger, math::get_multinomial},
};

#[builder]
pub fn resample(
    input: String,
    output: String,
    n_resamples: usize,
    seed: Option<u64>,
    samples: Vec<String>,
) -> Result<()> {
    let dataframe = load_dataframe(input.into())?;
    let count_matrix = dataframe
        .select(&samples)?
        .to_ndarray::<Float64Type>(IndexOrder::Fortran)?;

    let logger = Logger::new_silent();
    let normed_matrix = normalize_counts(&count_matrix, &Normalization::default(), &logger);
    let sgrna_counts = normed_matrix
        .mean_axis(Axis(1))
        .expect("Could not generate mean of sgRNAs across normed samples");
    let total = sgrna_counts.sum();
    let probabilities = sgrna_counts.mapv(|x| x / total);

    let mut resamples = vec![];
    let mut rng = if let Some(seed) = seed {
        ChaCha8Rng::seed_from_u64(seed)
    } else {
        ChaCha8Rng::from_entropy()
    };
    for idx in 0..n_resamples {
        let resample = get_multinomial(
            probabilities.as_slice().unwrap(),
            total.round() as u64,
            &mut rng,
        )?;
        let series = Series::new(&format!("resample_{idx}"), resample.to_vec());
        resamples.push(series);
    }

    let mut full_data = dataframe.hstack(&resamples)?;
    let writer = File::create(output).map(BufWriter::new)?;
    CsvWriter::new(writer)
        .with_separator(b'\t')
        .include_header(true)
        .with_quote_style(QuoteStyle::Never)
        .with_float_scientific(Some(true))
        .finish(&mut full_data)?;

    Ok(())
}
