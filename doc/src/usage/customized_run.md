# Command Arguments

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
crispr_screen --help
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

