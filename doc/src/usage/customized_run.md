# Command Arguments

## Subcommands

`crispr_screen` has two subcommands:
  1. `test`
  2. `agg`

`test` is used to perform the sgRNA-level differential abundance tests and then aggregate the results to the gene-level.
`test` by default will perform a gene-level aggregation as well, but can be skipped with the `--skip-agg` flag.

`agg` is used to just perform the gene-level aggregation on a precalculated differential abundance matrix.
Take a look at the `results.sgrna.tsv` file to see the expected file format required. Column names can
be provided as well - details can be found by running `crispr_screen agg --help`

## Arguments

### Required

The three required arguments are as follows:

| Argument | Description |
|-|-|
| **input** | Filepath of the input count matrix |
| **controls** | Labels for the control samples (can take multiple space-separated values) |
| **treatments** | Labels for the treatment samples (can take multiple space-separated values) |


### Optional Arguments

These arguments are optional and may change the configuration of the analysis.

For further details on them please run:

```bash
crispr_screen test --help
```

| Argument | Description |
|-|-|
| **output** | Prefix of the output sgRNA and gene result dataframes |
| **norm** | Normalization method to use |
| **agg** | Gene aggregation method to use |
| **correction** | Multiple hypothesis correction to use |
| **model-choice** | Which least squares model to fit |
| **alpha** | The alpha threshold parameter for aRRA algorithm |
| **permutations** | The number of permutations to perform in aRRA algorithm |
| **no-adjust-alpha** | Use flag to have fixed alpha, otherwise an empirical one will be calculated from provided alpha. |
| **ntc-token** | The token string to search for non-targeting controls (if INC) |

