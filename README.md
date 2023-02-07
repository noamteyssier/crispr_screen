# CRISPR Screen

[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](./LICENSE.md)
![actions status](https://github.com/noamteyssier/crispr_screen/workflows/CI/badge.svg)
[![codecov](https://codecov.io/gh/noamteyssier/crispr_screen/branch/main/graph/badge.svg?token=9ALCE60W2T)](https://codecov.io/gh/noamteyssier/crispr_screen)

This tool is recreation of the MAGeCK Differential Expression algorithm described
in the original paper with options included.
It also extends the algorithm and provides an option for the MAGeCK-INC pipeline
which performs a Mann-Whitney U-Test on the Non-Targeting Controls for gene score
aggregation as opposed to the traditional alpha-RRA algorithm.
The goal of this project is to have an efficient and clean implementation of a
CRISPR-screening platform that is as easy to use as possible without sacrificing
efficiency and interpretability.

## Installation

### Installing Rust

If you do not have the rust toolchain installed you can do so with the following
command

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### Installing CRISPR Screen

```bash
cargo install crispr_screen
```

## Usage

This tool is meant to run from the commandline and expects at minimum 3 arguments:

1. an sgRNA-Gene count table
2. the control label(s)
3. the treatment label(s).

The following list of examples is not exhaustive and if you are confused by any
parameters I recommend running the help menu:

```bash
crispr_screen --help
```

### Basic Run with a Single Control and a Single Treatment

Here is the command layout for a basic run with only these arguments:

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label>
```

### Basic Run with a Single Control and Multiple Treatments

If you have multiple treatment and/or control labels you can provide multiple
arguments to the same flag.
In this case I am providing two treatment labels for the same `-t` flag.

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label_a> <treatment_label_b>
```

### Basic Run with a Specified Output Prefix

By default this will write the sgRNA-level and Gene-level results to
`results.sgrna_results.tab` and `results.gene_results.tab` respectively,
but you can specificy your own output prefix (i.e. replacing `results`)
with the `-o` flag.

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label> -o my_prefix
```

### Alternative Normalization Strategies

Currently there are two normalization strategies provided by this tool -
`MedianRatio` and `Total`.
You can specify which one you'd like to use with the `-n` flag.
The two accepted options are `median-ratio` and `total`.

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label> -n total
```

### Alternative Gene-Aggregation Strategies

Currently there are two gene aggregation strategies provided by this tool -
`AlphaRRA` and `INC`.
By default the gene aggregation strategy run is `AlphaRRA`.
However, you can specify which one you'd like to use with the `-g` flag.
The two accepted options are `rra` and `inc`.

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label> -g inc
```
