# CRISPR Screen

This tool is faithful recreation of the MAGeCK Differential Expression algorithm described in the original paper. 
It also extends the algorithm and provides an option for the MAGeCK-INC pipeline which performs a Mann-Whitney U-Test on the Non-Targeting Controls for gene score aggregation as opposed to the traditional alpha-RRA algorithm.
This is a work in progress and new features will be added over the coming weeks. 
The goal of this project is to have an efficient and clean implementation of a CRISPR-screening platform that is as easy to use as possible without sacrificing efficiency and interpretability. 

## Installation

### Installing Rust
If you do not have the rust toolchain installed you can do so with the following command

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

If you do have rust installed I highly recommend updating to the newest toolchain to avoid any potential versioning issues

```bash
rustup update
```

### Installing CRISPR Screen
When this tool is stable I will upload to [crates.io](https://crates.io) for easy installation but until then you can install from the github repo.

```bash
git clone github.com/noamteyssier/crispr_screen
cd crispr_screen
cargo install --path .
```

## Usage
This tool is meant to run from the commandline and expects at minimum 3 arguments: an sgRNA-Gene count table, the control label(s), and the treatment label(s). 
The following list of examples is not exhaustive and if you are confused by any parameters I recommend running the help menu:
```bash
crispr_screen --help
```

### Basic Run with a Single Control and a Single Treatment
Here is the command layout for a basic run with only these arguments: 

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label>
```

### Basic Run with a Single Control and Multiple Treatments
If you have multiple treatment and/or control labels you can provide multiple arguments to the same flag. 
In this case I am providing two treatment labels for the same `-t` flag.

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label_a> <treatment_label_b>
```

### Basic Run with a Specified Output Prefix
By default this will write the sgRNA-level and Gene-level results to `results.sgrna_results.tab` and `results.gene_results.tab` respectively, but you can specificy your own output prefix (i.e. replacing `results`) with the `-o` flag.

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label> -o my_prefix
```


### Alternative Normalization Strategies.
Currently there are two normalization strategies provided by this tool - `MedianRatio` and `Total`.
You can specify which one you'd like to use with the `-n` flag.
The two accepted options are `median` and `total`.

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label> -n total
```

### Alternative Gene-Aggregation Strategies.
Currently there are two gene aggregation strategies provided by this tool - `AlphaRRA` and `INC`.
By default the gene aggregation strategy run is `AlphaRRA`.
However, you can specify which one you'd like to use with the `-g` flag.
The two accepted options are `rra` and `inc`.

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label> -g inc
```

### RRA Specific Options
There are two `AlphaRRA` specific parameterizations. 
The first is the `alpha` parameter which is used to paramterize the p-value threshold to consider when doing the permutation tests. 
This is set to `0.1` by default, but you can override that parameterization with the `-a` flag.
The second is the `permutations` paramater which is then multiplied by the number of genes to determines the number of total permutations to test.
This is set to `100` by default (representing `100 * N_GENES` total permutations), but you can override that with the `-p` flag.

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label> -a 0.05 -p 500
```

### INC Specific Options
The `INC` algorithm requires knowing which of your sgRNAs are non-targeting controls. 
The way this is done is by searching the associated gene for the control-token which by default is set to `non-targeting`.
All genes which contain the control-token will be considering non-targeting controls and used as the null distribution when performing Mann-Whitney U Tests.
If your control-token is different however or you'd like to specify a novel one you can do so with the `-n` flag.

```bash
crispr_screen -i <count_table> -c <control_label> -t <treatment_label> -g inc -n NTC
```

