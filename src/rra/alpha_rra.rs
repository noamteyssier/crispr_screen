use ndarray::Array1;
use super::{encode_index, normed_ranks, group_sizes};

pub fn alpha_rra(
    pvalues: &Array1<f64>,
    genes: &Vec<String>) 
{
    let encode = encode_index(genes);
    let nranks = normed_ranks(pvalues);
    let sizes = group_sizes(&encode);
}

