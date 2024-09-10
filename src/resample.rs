use std::{fs::File, io::BufWriter};

use anyhow::Result;
use bon::builder;
use ndarray::prelude::*;
use ndarray_rand::rand::SeedableRng;
use polars::prelude::*;
use rand_chacha::ChaCha8Rng;
use rand_distr::{Dirichlet, DirichletError, Distribution};

use crate::{
    io::{load_dataframe, to_ndarray},
    norm::{normalize_counts, Normalization},
    utils::{logging::Logger, math::get_multinomial},
};

/// Builds a Dirichlet distribution from a given array of alpha values
///
/// Adds 1 to each value in the array to ensure that the distribution is valid
///
/// ```text
/// α ~ Dirichlet(x + 1)
/// ```
fn build_dirichlet(alpha: &Array1<f64>) -> Result<Dirichlet<f64>, DirichletError> {
    let alpha: Vec<f64> = alpha.iter().map(|x| x + 1.0).collect();
    Dirichlet::new(&alpha)
}

/// Resamples the counts of a given Dirichlet distribution
///
/// This function will sample from the Dirichlet distribution and then resample the counts
/// based on the total sum of the counts using a multinomial distribution
///
/// ```text
/// p ~ Dirichlet(alpha)
/// resample ~ Multinomial(p, total)
/// ```
fn resample_counts(
    dirichlet: &Dirichlet<f64>,
    total: f64,
    rng: &mut ChaCha8Rng,
) -> Result<Vec<u64>> {
    let probabilities = dirichlet.sample(rng);
    let resample = get_multinomial(&probabilities, total.round() as u64, rng)?;
    Ok(resample)
}

/// Loops over the resampling process to generate multiple resamples
/// of the counts from a Dirichlet-Multinomial distribution
fn loop_resample(
    dirichlet: &Dirichlet<f64>,
    total: f64,
    rng: &mut ChaCha8Rng,
    n_resamples: usize,
) -> Result<Vec<Series>> {
    let mut resamples = vec![];
    for idx in 0..n_resamples {
        let resample = resample_counts(dirichlet, total, rng)?;
        let series = Series::new(&format!("resample_{idx}"), resample);
        resamples.push(series);
    }
    Ok(resamples)
}

/// Builds a ChaCha8Rng from a given seed or entropy
fn build_rng(seed: Option<u64>) -> ChaCha8Rng {
    if let Some(seed) = seed {
        ChaCha8Rng::seed_from_u64(seed)
    } else {
        ChaCha8Rng::from_entropy()
    }
}

#[builder]
pub fn resample(
    input: String,
    output: String,
    n_resamples: usize,
    seed: Option<u64>,
    samples: Vec<String>,
) -> Result<()> {
    let dataframe = load_dataframe(input.into())?;
    let count_matrix = to_ndarray(&dataframe, &samples)?;

    let logger = Logger::new_silent();
    let normed_matrix = normalize_counts(&count_matrix, &Normalization::default(), &logger);
    let sgrna_counts = normed_matrix
        .mean_axis(Axis(1))
        .expect("Could not generate mean of sgRNAs across normed samples");
    let total = sgrna_counts.sum();

    let mut rng = build_rng(seed);
    let dirichlet = build_dirichlet(&sgrna_counts)?;
    let resamples = loop_resample(&dirichlet, total, &mut rng, n_resamples)?;

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
