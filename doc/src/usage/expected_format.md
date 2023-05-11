# Expected Formats

## Inputs

### Count Table

In the above example, the `count_table.tsv` is a **tab-separated file** of the
following form:

| sgrna | gene | low_1 | high_1 | low_2 | high_2 |
|-------|------|-------|--------|-------|--------|
| sgrna.0 | gene.0 | 1 | 4 | 2 | 6 |
| sgrna.1 | gene.0 | 2 | 6 | 3 | 4 |
| ...     | ...    | ... | ... | ... | ... |
| sgrna.n | gene.m | 4 | 12 | 5 | 20 |

> **_Note_**:
>
> The sample names provided **do not** need to be in the same order you
> they appear in the file.
>
> However, the first two columns **do need** to be an sgRNA and gene column
> respectively (though they can be named whatever you like.)

If you have extra columns in the table you don't want to analyze, just provide
the names of the columns you do want to analyze.

The above command will still work as expected for the following table:

| sgrna | gene | low_1 | high_1 | low_2 | high_2 | ext_1 |
|-------|------|-------|--------|-------|--------|-------|
| sgrna.0 | gene.0 | 1 | 4 | 2 | 6 | 20 |
| sgrna.1 | gene.0 | 2 | 6 | 3 | 4 | 30 |
| ...     | ...    | ... | ... | ... | ... | ... |
| sgrna.n | gene.m | 4 | 12 | 5 | 20 | 5 |

## Outputs

### sgRNA Results

The sgRNA results dataframe (written to `<args.output>.sgrna_results.tsv`) is a
table whose columns are of the following form:

| Column | Description |
|--------|-------------|
| **sgrna** | The sgRNA name provided in the first column of the `count_table`. |
| **gene** | The gene name provided in the second column of the `count_table`. |
| **base** | The normalized/aggregated count of the sgRNA across all samples. |
| **control** | The normalized/aggregated count of the sgRNA for the controls. |
| **treatment** | The normalized/aggregated count of the sgRNA for the treatment. |
| **adj_var** | The adjusted variance for the sgRNA determined by the least squares fit. |
| **fold_change** | The fold change of the treatment over the controls. |
| **log2_fold_change** | The log2 fold change of the treatment over the controls. |
| **pvalue_low** | The p-value for an depletion of the sgRNA. |
| **pvalue_high** | The p-value for an enrichment of the sgRNA. |
| **pvalue_twosided** | The two-sided p-value of an enrichment or depletion of the sgRNA. |
| **fdr** | The adjusted false discovery rate of the sgRNA. |

### Gene Results

The gene results dataframe (written to `<args.output>.gene_results.tsv`) is a
table whose columns are of the following form:

| Column | Description |
|--------|-------------|
| **gene** | The gene name provided in the second column of the `count_table`. |
| **fold_change** | The aggregated fold change of the treatment from the controls. |
| **log_fold_change** | The log2 aggregated fold change of the treatment from the controls. |
| **score_low** | The minimum p-value observed in the RRA or the U-score observed in the INC for the gene being depleted. |
| **pvalue_low** | The aggregated p-value for a depletion of the gene. |
| **fdr_low** | The false discovery rate for a depletion of the gene. |
| **score_high** | The minimum p-value observed in the RRA of the U-score observed in the INC for the gene being enriched. |
| **pvalue_high** | The aggregated p-value for an enrichment of the gene. |
| **fdr_high** | The false discovery rate for an enrichment of the gene. |
| **pvalue** | The minimum pvalue observed with either test. |
| **fdr** | The minimum false discovery rate observed with either test. |
| **phenotype_score** | The -log10 FDR multiplied by the log2 fold change of the gene. |

> Note: If you ran `crispr_screen` with `INC`
>
> The FDR in this case will be the **empirical false discovery rate** given the non-targeting controls.
> As a result it will **not** be strictly monotonic with the p-values.
> I recommend working with the pvalues in this case and paying attention to the calculated thresholds
> of the FDR given in the workflow log.

### Hit Results

The hits dataframe (written to `<args.output>.hits.tsv`) is a
table whose columns are of the following form:

| Column | Description |
|--------|-------------|
| **gene** | The gene name provided in the second column of the `count_table`. |
| **log2fc** | The log2 aggregated fold change of the treatment from the controls. |
| **pvalue** | The minimum p-value observed in the aggregation test (minimum of both sides). |
| **phenotype_score** | The product of the `log2fc` and the `-log10(pvalue)`. |
| **fdr** | The calculated false discovery rate (only shown if running $\alpha$-RRA). |
