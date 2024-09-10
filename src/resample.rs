use anyhow::{bail, Result};
use bon::builder;
use ndarray::prelude::*;
use ndarray_rand::rand::SeedableRng;
use polars::prelude::*;
use rand_chacha::ChaCha8Rng;
use rand_distr::{Dirichlet, DirichletError, Distribution};

use crate::{
    io::{build_regex_set, load_dataframe, match_headers_from_regex_set, to_ndarray, write_tsv},
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

fn calculate_depth(
    dataframe: &DataFrame,
    sgrna_counts: &Array1<f64>,
    depth: Option<usize>,
    depth_samples: Option<Vec<String>>,
) -> Result<f64> {
    if depth.is_some() & depth_samples.is_some() {
        bail!("Cannot specify both depth and depth_samples at the same time")
    }
    if let Some(depth) = depth {
        Ok(depth as f64)
    } else if let Some(depth_samples) = depth_samples {
        let sample_regex = build_regex_set(&depth_samples)?;
        let sample_labels = match_headers_from_regex_set(dataframe, &sample_regex)?;
        let count_matrix = to_ndarray(dataframe, &sample_labels)?;
        let sample_sums = count_matrix.sum_axis(Axis(0));
        let mean_depth = sample_sums.mean().expect("Could not calculate mean depth");
        Ok(mean_depth)
    } else {
        Ok(sgrna_counts.iter().sum())
    }
}

#[builder]
pub fn resample(
    input: String,
    path: Option<String>,
    n_resamples: usize,
    depth: Option<usize>,
    depth_samples: Option<Vec<String>>,
    seed: Option<u64>,
    samples: Vec<String>,
    quiet: bool,
) -> Result<()> {
    let logger = Logger::from_quiet(quiet);
    logger.start_resampling();

    let dataframe = load_dataframe(input.into())?;

    let sample_regex = build_regex_set(&samples)?;
    let sample_labels = match_headers_from_regex_set(&dataframe, &sample_regex)?;
    logger.sampled_names(&sample_labels);
    logger.number_of_resamples(n_resamples);

    let count_matrix = to_ndarray(&dataframe, &sample_labels)?;
    let normed_matrix = normalize_counts(&count_matrix, &Normalization::default(), &logger);
    let sgrna_counts = normed_matrix
        .mean_axis(Axis(1))
        .expect("Could not generate mean of sgRNAs across normed samples");
    logger.num_sgrnas(sgrna_counts.as_slice().expect("Could not create slice"));

    let total = calculate_depth(&dataframe, &sgrna_counts, depth, depth_samples)?;
    logger.resampling_depth(total as usize);

    let mut rng = build_rng(seed);
    logger.describe_seed(seed);
    let dirichlet = build_dirichlet(&sgrna_counts)?;
    let resamples = loop_resample(&dirichlet, total, &mut rng, n_resamples)?;

    let mut full_data = dataframe.hstack(&resamples)?;
    write_tsv(&mut full_data, path)?;

    Ok(())
}
