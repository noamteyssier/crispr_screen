use std::{
    fmt::{Display, Formatter},
    fs::File,
    io::{BufWriter, Write},
};

use crate::{
    aggregation::{AggregationResult, GeneAggregation},
    utils::config::Configuration,
};
use anyhow::Result;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
enum Method {
    AlphaRRA,
    IncProd,
    IncPvalue,
    GeoPAGG,
}
impl Display for Method {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Method::AlphaRRA => write!(f, "rra"),
            Method::IncProd => write!(f, "inc-product"),
            Method::IncPvalue => write!(f, "inc-pvalue"),
            Method::GeoPAGG => write!(f, "geopagg"),
        }
    }
}

pub struct Screenviz {
    method: Method,
    gene: String,
    x: String,
    y: String,
    z: String,
    threshold: Option<f64>,
    threshold_low: Option<f64>,
    threshold_high: Option<f64>,
    ntc_token: Option<String>,
}
impl Screenviz {
    pub fn new(results: &AggregationResult, config: &Configuration) -> Self {
        let method = match config.aggregation() {
            GeneAggregation::AlpaRRA {
                alpha: _,
                npermutations: _,
                adjust_alpha: _,
                fdr: _,
            } => Method::AlphaRRA,
            GeneAggregation::Inc {
                token: _,
                fdr: _,
                group_size: _,
                n_draws: _,
                use_product,
            } => {
                if *use_product {
                    Method::IncProd
                } else {
                    Method::IncPvalue
                }
            }
            GeneAggregation::GeoPAGG {
                token: _,
                weight_config: _,
                fdr: _,
            } => Method::GeoPAGG,
        };
        let gene = "gene".to_string();
        let x = "log_fold_change".to_string();
        let y = "pvalue".to_string();
        let z = match config.aggregation() {
            GeneAggregation::AlpaRRA {
                alpha: _,
                npermutations: _,
                adjust_alpha: _,
                fdr: _,
            } => "fdr".to_string(),
            GeneAggregation::Inc {
                token: _,
                fdr: _,
                group_size: _,
                n_draws: _,
                use_product,
            } => {
                if *use_product {
                    "phenotype_score".to_string()
                } else {
                    "pvalue".to_string()
                }
            }
            GeneAggregation::GeoPAGG {
                token: _,
                weight_config: _,
                fdr: _,
            } => "fdr".to_string(),
        };

        let (threshold, threshold_low, threshold_high, ntc_token) = match config.aggregation() {
            GeneAggregation::AlpaRRA {
                alpha: _,
                npermutations: _,
                adjust_alpha: _,
                fdr,
            } => (Some(*fdr), None, None, None),
            GeneAggregation::Inc {
                token: _,
                fdr: _,
                group_size: _,
                n_draws: _,
                use_product: _,
            } => (
                None,
                results.threshold_low(),
                results.threshold_high(),
                Some(String::from("pseudogene")),
            ),
            GeneAggregation::GeoPAGG {
                token: _,
                weight_config: _,
                fdr,
            } => (Some(*fdr), None, None, Some(String::from("amalgam"))),
        };

        Self {
            method,
            gene,
            x,
            y,
            z,
            threshold,
            threshold_low,
            threshold_high,
            ntc_token,
        }
    }

    pub fn write(&self, prefix: &str) -> Result<()> {
        let mut writer = File::create(format!("{prefix}.screenviz.yaml")).map(BufWriter::new)?;
        match self.method {
            Method::AlphaRRA => {
                writeln!(writer, "method: {}", self.method)?;
                writeln!(writer, "gene: {}", self.gene)?;
                writeln!(writer, "x: {}", self.x)?;
                writeln!(writer, "y: {}", self.y)?;
                writeln!(writer, "z: {}", self.z)?;
                writeln!(writer, "threshold: {}", self.threshold.unwrap())?;
            }
            Method::GeoPAGG => {
                writeln!(writer, "method: {}", self.method)?;
                writeln!(writer, "gene: {}", self.gene)?;
                writeln!(writer, "x: {}", self.x)?;
                writeln!(writer, "y: {}", self.y)?;
                writeln!(writer, "z: {}", self.z)?;
                writeln!(writer, "threshold: {}", self.threshold.unwrap())?;
                writeln!(writer, "ntc_token: {}", self.ntc_token.as_ref().unwrap())?;
            }
            _ => {
                writeln!(writer, "method: {}", self.method)?;
                writeln!(writer, "gene: {}", self.gene)?;
                writeln!(writer, "x: {}", self.x)?;
                writeln!(writer, "y: {}", self.y)?;
                writeln!(writer, "z: {}", self.z)?;
                writeln!(writer, "threshold_low: {}", self.threshold_low.unwrap())?;
                writeln!(writer, "threshold_high: {}", self.threshold_high.unwrap())?;
                writeln!(writer, "ntc_token: {}", self.ntc_token.as_ref().unwrap())?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod testing {
    use super::*;
    use ndarray::Array1;

    fn build_config_rra<'a>() -> Configuration<'a> {
        let aggregation = GeneAggregation::AlpaRRA {
            alpha: 0.05,
            npermutations: 1000,
            adjust_alpha: true,
            fdr: 0.05,
        };
        let base_mean = 0.0;
        let prefix = "results";
        Configuration::builder()
            .aggregation(aggregation)
            .min_base_mean(base_mean)
            .prefix(prefix)
            .build()
    }

    fn build_config_prod<'a>() -> Configuration<'a> {
        let aggregation = GeneAggregation::Inc {
            token: "non-targeting",
            fdr: 0.05,
            group_size: 5,
            use_product: true,
            n_draws: 100,
        };
        let prefix = "results";
        Configuration::builder()
            .aggregation(aggregation)
            .prefix(prefix)
            .build()
    }

    fn build_config_pvalue<'a>() -> Configuration<'a> {
        let aggregation = GeneAggregation::Inc {
            token: "non-targeting",
            fdr: 0.05,
            group_size: 5,
            use_product: false,
            n_draws: 100,
        };
        let prefix = "results";
        Configuration::builder()
            .aggregation(aggregation)
            .prefix(prefix)
            .build()
    }

    fn build_results_rra() -> AggregationResult {
        let genes = vec!["gene1".to_string(), "gene2".to_string()];
        let gene_fc = Array1::from(vec![1.0, 2.0]);
        let pvalues_low = Array1::from(vec![0.1, 0.2]);
        let pvalues_high = Array1::from(vec![0.3, 0.4]);
        let fdr_low = Array1::from(vec![0.5, 0.6]);
        let fdr_high = Array1::from(vec![0.7, 0.2]);
        let aggregation_score_low = Array1::from(vec![0.5, 0.6]);
        let aggregation_score_high = Array1::from(vec![0.7, 0.8]);
        let threshold_low = None;
        let threshold_high = None;
        AggregationResult::new(
            genes,
            gene_fc,
            pvalues_low,
            pvalues_high,
            fdr_low,
            fdr_high,
            aggregation_score_low,
            aggregation_score_high,
            threshold_low,
            threshold_high,
        )
    }

    fn build_results_incprod() -> AggregationResult {
        let genes = vec!["gene1".to_string(), "gene2".to_string()];
        let gene_fc = Array1::from(vec![1.0, 2.0]);
        let pvalues_low = Array1::from(vec![0.1, 0.2]);
        let pvalues_high = Array1::from(vec![0.3, 0.4]);
        let fdr_low = Array1::from(vec![0.5, 0.6]);
        let fdr_high = Array1::from(vec![0.7, 0.2]);
        let aggregation_score_low = Array1::from(vec![0.5, 0.6]);
        let aggregation_score_high = Array1::from(vec![0.7, 0.8]);
        let threshold_low = Some(-0.5);
        let threshold_high = Some(0.5);
        AggregationResult::new(
            genes,
            gene_fc,
            pvalues_low,
            pvalues_high,
            fdr_low,
            fdr_high,
            aggregation_score_low,
            aggregation_score_high,
            threshold_low,
            threshold_high,
        )
    }

    fn build_results_incpvalue() -> AggregationResult {
        let genes = vec!["gene1".to_string(), "gene2".to_string()];
        let gene_fc = Array1::from(vec![1.0, 2.0]);
        let pvalues_low = Array1::from(vec![0.1, 0.2]);
        let pvalues_high = Array1::from(vec![0.3, 0.4]);
        let fdr_low = Array1::from(vec![0.5, 0.6]);
        let fdr_high = Array1::from(vec![0.7, 0.2]);
        let aggregation_score_low = Array1::from(vec![0.5, 0.6]);
        let aggregation_score_high = Array1::from(vec![0.7, 0.8]);
        let threshold_low = Some(0.005);
        let threshold_high = Some(0.004);
        AggregationResult::new(
            genes,
            gene_fc,
            pvalues_low,
            pvalues_high,
            fdr_low,
            fdr_high,
            aggregation_score_low,
            aggregation_score_high,
            threshold_low,
            threshold_high,
        )
    }

    #[test]
    fn test_screenviz_incprod() {
        let config = build_config_prod();
        let results = build_results_incprod();
        let screenviz = Screenviz::new(&results, &config);
        assert_eq!(screenviz.method, Method::IncProd);
        assert_eq!(screenviz.gene, "gene");
        assert_eq!(screenviz.x, "log_fold_change");
        assert_eq!(screenviz.y, "pvalue");
        assert_eq!(screenviz.z, "phenotype_score");
        assert_eq!(screenviz.threshold, None);
        assert_eq!(screenviz.threshold_low, Some(-0.5));
        assert_eq!(screenviz.threshold_high, Some(0.5));
        assert_eq!(screenviz.ntc_token, Some("pseudogene".to_string()));
    }

    #[test]
    fn test_screenviz_incpvalue() {
        let config = build_config_pvalue();
        let results = build_results_incpvalue();
        let screenviz = Screenviz::new(&results, &config);
        assert_eq!(screenviz.method, Method::IncPvalue);
        assert_eq!(screenviz.gene, "gene");
        assert_eq!(screenviz.x, "log_fold_change");
        assert_eq!(screenviz.y, "pvalue");
        assert_eq!(screenviz.z, "pvalue");
        assert_eq!(screenviz.threshold, None);
        assert_eq!(screenviz.threshold_low, Some(0.005));
        assert_eq!(screenviz.threshold_high, Some(0.004));
        assert_eq!(screenviz.ntc_token, Some("pseudogene".to_string()));
    }

    #[test]
    fn test_screenviz_rra() {
        let config = build_config_rra();
        let results = build_results_rra();
        let screenviz = Screenviz::new(&results, &config);
        assert_eq!(screenviz.method, Method::AlphaRRA);
        assert_eq!(screenviz.gene, "gene");
        assert_eq!(screenviz.x, "log_fold_change");
        assert_eq!(screenviz.y, "pvalue");
        assert_eq!(screenviz.z, "fdr");
        assert_eq!(screenviz.threshold, Some(0.05));
        assert_eq!(screenviz.threshold_low, None);
        assert_eq!(screenviz.threshold_high, None);
        assert_eq!(screenviz.ntc_token, None);
    }
}
