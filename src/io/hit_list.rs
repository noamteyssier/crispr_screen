use std::{
    fs::File,
    io::{BufWriter, Write},
};

use crate::{
    aggregation::{AggregationResult, GeneAggregation},
    utils::config::Configuration,
};
use anyhow::Result;
use hashbrown::HashSet;
use ndarray::{Array1, Axis};

#[derive(Clone, Copy)]
pub enum MethodEnum {
    /// Threshold set by product score
    IncProduct {
        threshold_low: f64,
        threshold_high: f64,
    },

    /// Threshold set by pvalue
    IncPvalue {
        threshold_low: f64,
        threshold_high: f64,
    },

    /// Threshold set by fdr
    RRA { fdr: f64 },
}
impl MethodEnum {
    pub fn new(result: &AggregationResult, config: &Configuration) -> Self {
        match config.aggregation() {
            GeneAggregation::Inc {
                token: _,
                fdr: _,
                group_size: _,
                n_draws: _,
                use_product,
            } => {
                if *use_product {
                    Self::IncProduct {
                        threshold_low: result.threshold_low().unwrap(),
                        threshold_high: result.threshold_high().unwrap(),
                    }
                } else {
                    Self::IncPvalue {
                        threshold_low: result.threshold_low().unwrap(),
                        threshold_high: result.threshold_high().unwrap(),
                    }
                }
            }
            GeneAggregation::AlpaRRA {
                alpha: _,
                npermutations: _,
                adjust_alpha: _,
                fdr,
            } => Self::RRA { fdr: *fdr },
        }
    }
}

#[derive(Clone, Copy)]
enum Direction {
    Less,
    Greater,
}

pub struct HitList {
    gene: Vec<String>,
    log2fc: Array1<f64>,
    pvalues: Array1<f64>,
    fdr: Option<Array1<f64>>,
    phenotype: Array1<f64>,
    size: usize,
    method: MethodEnum,
}

impl HitList {
    pub fn new(result: &AggregationResult, config: &Configuration) -> Self {
        let method = MethodEnum::new(result, config);
        let mask = Self::generate_mask(method, result);
        let gene = Self::select_genes(result, &mask);
        let log2fc = result.gene_log2_fc().select(Axis(0), &mask);
        let pvalues = result.pvalue().select(Axis(0), &mask);
        let phenotype = result.phenotype_score().select(Axis(0), &mask);
        let fdr = match method {
            MethodEnum::RRA { fdr: _ } => Some(result.fdr().select(Axis(0), &mask)),
            _ => None,
        };
        Self {
            size: gene.len(),
            gene,
            log2fc,
            pvalues,
            phenotype,
            fdr,
            method,
        }
    }

    fn generate_mask(method: MethodEnum, result: &AggregationResult) -> Vec<usize> {
        match method {
            MethodEnum::RRA { fdr } => Self::index_threshold(result.fdr(), fdr, Direction::Less),
            MethodEnum::IncProduct {
                threshold_low,
                threshold_high,
            } => {
                let mask_low =
                    Self::index_threshold(result.phenotype_score(), threshold_low, Direction::Less);
                let mask_high = Self::index_threshold(
                    result.phenotype_score(),
                    threshold_high,
                    Direction::Greater,
                );
                let mask: HashSet<usize> =
                    HashSet::from_iter(mask_low.iter().chain(mask_high.iter()).cloned());
                mask.into_iter().collect()
            }
            MethodEnum::IncPvalue {
                threshold_low,
                threshold_high,
            } => {
                let mask_low =
                    Self::index_threshold(result.pvalues_low(), threshold_low, Direction::Less);
                let mask_high =
                    Self::index_threshold(result.pvalues_high(), threshold_high, Direction::Less);
                let mask: HashSet<usize> =
                    HashSet::from_iter(mask_low.iter().chain(mask_high.iter()).cloned());
                mask.into_iter().collect()
            }
        }
    }

    fn index_threshold(array: &Array1<f64>, threshold: f64, lt: Direction) -> Vec<usize> {
        match lt {
            Direction::Less => array
                .iter()
                .enumerate()
                .filter(|(_, &x)| x <= threshold)
                .map(|(i, _)| i)
                .collect(),
            Direction::Greater => array
                .iter()
                .enumerate()
                .filter(|(_, &x)| x >= threshold)
                .map(|(i, _)| i)
                .collect(),
        }
    }

    fn select_genes(result: &AggregationResult, mask: &[usize]) -> Vec<String> {
        let mut genes = Vec::with_capacity(mask.len());
        for i in mask.iter() {
            genes.push(result.genes()[*i].clone());
        }
        genes
    }

    pub fn write(&self, prefix: &str) -> Result<()> {
        let mut writer = File::create(format!("{}.hits.tsv", prefix)).map(BufWriter::new)?;
        match self.method {
            MethodEnum::RRA { fdr: _ } => {
                writeln!(writer, "gene\tlog2fc\tpvalue\tphenotype_score\tfdr")?;

                for idx in 0..self.size {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}\t{}",
                        self.gene[idx],
                        self.log2fc[idx],
                        self.pvalues[idx],
                        self.phenotype[idx],
                        self.fdr.as_ref().unwrap()[idx]
                    )?;
                }
            }
            _ => {
                writeln!(writer, "gene\tlog2fc\tpvalue\tphenotype_score")?;

                for idx in 0..self.size {
                    writeln!(
                        writer,
                        "{}\t{}\t{}\t{}",
                        self.gene[idx], self.log2fc[idx], self.pvalues[idx], self.phenotype[idx],
                    )?;
                }
            }
        };
        Ok(())
    }

    pub fn num_total(&self) -> usize {
        self.size
    }

    pub fn num_enrichments(&self) -> usize {
        self.log2fc.iter().filter(|&&x| x > 0.0).count()
    }

    pub fn num_depletions(&self) -> usize {
        self.log2fc.iter().filter(|&&x| x < 0.0).count()
    }
}
