use hashbrown::HashMap;
use ndarray::Array1;
use super::{
    normed_ranks, group_sizes, 
    permutations::run_permutations, filter_alpha, 
    robust_rank::robust_rank_aggregation, 
    utils::empirical_cdf};
use crate::{aggregation::utils::{encode_index, select_ranks}, utils::logging::Logger};

/// Calculates an empirical p-value of the robust rank aggregation for the current gene set with
/// respect to random permutations of that size
fn gene_rra(
    current_idx: usize,
    encodings: &[usize],
    nranks: &Array1<f64>,
    permutation_vectors: &HashMap<usize, Array1<f64>>,
    alpha: f64) -> (f64, f64)
{
    let gene_ranks = select_ranks(current_idx, encodings, nranks);
    let filtered = filter_alpha(&gene_ranks, alpha);
    let score = robust_rank_aggregation(&filtered, gene_ranks.len());
    (
        score,
        empirical_cdf(score, permutation_vectors.get(&gene_ranks.len()).expect("Unexpected missing key"))
    )
}

/// Performs the alpha-RRA algorithm
pub fn alpha_rra(
    pvalues: &Array1<f64>,
    genes: &Vec<String>,
    alpha: f64,
    npermutations: usize,
    logger: &Logger) -> (Vec<String>, Array1<f64>, Array1<f64>)
{
    let (encode_map, encode) = encode_index(genes);
    let n_genes = encode_map.len();
    let nranks = normed_ranks(pvalues);
    let sizes = group_sizes(&encode);
    logger.permutation_sizes(&sizes);

    // calculate rra scores for a vector of random samplings for each unique size
    let permutation_vectors = sizes
        .iter()
        .map(|unique_size| (*unique_size, run_permutations(nranks.len(), alpha, npermutations * n_genes, *unique_size)))
        .map(|(u, v)| (u, Array1::from_vec(v)))
        .collect::<HashMap<usize, Array1<f64>>>();

    // calculate empirical pvalues for each of the gene sets given the random nulls
    let (scores, pvalues): (Vec<f64>, Vec<f64>) = (0..*encode.iter().max().expect("Unexpected empty encoding"))
        .map(|curr| gene_rra(curr, &encode, &nranks, &permutation_vectors, alpha))
        .unzip();

    let names = (0..*encode.iter().max().expect("Unexpected empty encoding"))
        .map(|curr| encode_map.get(&curr).expect("Unexpected missing index").clone())
        .collect();

    (
        names, 
        Array1::from_vec(scores),
        Array1::from_vec(pvalues)
    )
}
