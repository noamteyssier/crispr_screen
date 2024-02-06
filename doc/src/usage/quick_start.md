# Quick Start

## Input Files

`crispr_screen` is meant to be run from the command-line and expects **at minimum**
three arguments:

1. an sgRNA-Gene count table (see [formats](./expected_format.md) for details)
2. the labels of the controls
3. the labels of the treatments

For additional options and descriptions see [arguments](./customized_run.md).

## Running [ `crispr_screen` ]

```bash
crispr_screen test \
  -i count_table.tsv \
  -c low_1 low_2 \
  -t high_1 high_2
```

## Output Files

This will create **three results files**:

1. `results.sgrna_results.tsv`
2. `results.gene_results.tsv`
3. `results.hits.tsv`

With differential expression statistics for each respectively
and a table for just the genes that are found to be significant
(for descriptions of files see [formats](./expected_format.md)).
