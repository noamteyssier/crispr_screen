# Methods

## Overview

The [ `crispr_screen` ] pipeline is really two analyses in one.

There is first the **sgRNA enrichment** analysis, and second the
**sgRNA aggregation** analysis.

The reasoning for this is that we want to be able to determine
which sgRNAs are statistically enriched/depleted in the treatments
given the controls from their abundances. However, because there
are generally such a high number of sgRNAs, there will be likely
be many sgRNA *hits* that cross spuriously either statistically
*or* biologically. If we were to just say that any gene with a
statistically significant sgRNA was significantly enriched or
depleted we will likely not end up with a set of genes which are
*robust* in our assay.

### sgRNA Enrichment

An overview of the sgRNA enrichment analysis is as follows:

- [Normalize the read counts](./methods/normalization.md)
- [Fit the dispersion model](./methods/dispersion.md)
- [Perform sgRNA-level enrichment](./methods/enrichment.md)
- [Adjust p-values for multiple hypothesis correction](./methods/correction.md)

### sgRNA Aggregation

An overview of the sgRNA aggregation procedure is as follows:

- Perform sgRNA aggregation (either `αRRA` or `INC`)
  - [αRRA](./methods/rra.md)
  - [INC](./methods/inc.md)
- [Adjust p-values for multiple hypothesis correction](./methods/correction.md)
