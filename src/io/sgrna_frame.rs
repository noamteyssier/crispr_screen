use crate::enrich::EnrichmentResult;
use anyhow::Result;
use ndarray::Array1;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

pub struct SgrnaFrame<'a> {
    sgrna_names: &'a [String],
    gene_names: &'a [String],
    adj_var: &'a Array1<f64>,
    base: &'a Array1<f64>,
    control: &'a Array1<f64>,
    treatment: &'a Array1<f64>,
    pvalue_low: &'a Array1<f64>,
    pvalue_high: &'a Array1<f64>,
    pvalue_twosided: &'a Array1<f64>,
    fdr: &'a Array1<f64>,
    fold_change: &'a Array1<f64>,
    log_fold_change: &'a Array1<f64>,
    product: &'a Array1<f64>,
    size: usize,
}
impl<'a> SgrnaFrame<'a> {
    pub fn new(
        sgrna_names: &'a [String],
        gene_names: &'a [String],
        adj_var: &'a Array1<f64>,
        sgrna_results: &'a EnrichmentResult,
    ) -> Self {
        Self {
            sgrna_names,
            gene_names,
            adj_var,
            base: sgrna_results.base_means(),
            control: sgrna_results.control_means(),
            treatment: sgrna_results.treatment_means(),
            pvalue_low: sgrna_results.pvalues_low(),
            pvalue_high: sgrna_results.pvalues_high(),
            pvalue_twosided: sgrna_results.pvalues_twosided(),
            fdr: sgrna_results.fdr(),
            fold_change: sgrna_results.fold_change(),
            log_fold_change: sgrna_results.log_fold_change(),
            size: sgrna_names.len(),
            product: sgrna_results.product(),
        }
    }

    pub fn write(&self, prefix: &str) -> Result<()> {
        let mut writer = File::create(format!("{prefix}.sgrna_results.tsv")).map(BufWriter::new)?;

        writeln!(
            writer,
            "sgrna\tgene\tbase\tcontrol\ttreatment\tadj_var\tfold_change\tlog2_fold_change\tpvalue_low\tpvalue_high\tpvalue_twosided\tfdr\tproduct",
        )?;

        for idx in 0..self.size {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}",
                self.sgrna_names[idx],
                self.gene_names[idx],
                self.base[idx],
                self.control[idx],
                self.treatment[idx],
                self.adj_var[idx],
                self.fold_change[idx],
                self.log_fold_change[idx],
                self.pvalue_low[idx],
                self.pvalue_high[idx],
                self.pvalue_twosided[idx],
                self.fdr[idx],
                self.product[idx],
            )?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod testing {
    use adjustp::Procedure;

    use crate::enrich::EnrichmentResult;

    #[test]
    fn test_sgrna_frame() {
        let sgrna_names = vec!["sgrna1".to_string(), "sgrna2".to_string()];
        let gene_names = vec!["gene1".to_string(), "gene2".to_string()];
        let adj_var = vec![0.1, 0.2].into_iter().collect();
        let control = vec![1.0, 2.0].into_iter().collect();
        let treatment = vec![3.0, 4.0].into_iter().collect();
        let pvalue_low = vec![0.1, 0.2].into_iter().collect();
        let pvalue_high = vec![0.3, 0.4].into_iter().collect();
        let correction = Procedure::BenjaminiHochberg;

        let enrichment =
            EnrichmentResult::new(pvalue_low, pvalue_high, control, treatment, correction);

        let sgrna_frame = super::SgrnaFrame::new(&sgrna_names, &gene_names, &adj_var, &enrichment);
        assert_eq!(sgrna_frame.size, 2);
        sgrna_frame.write("test").unwrap();
    }
}
