
## Description

Represents the methylation levels and the detected DMRs present in a defined region of the genome for a set of samples.

## Dependencies

``` r
# Gviz
BiocManager::install("Gviz")

# hg19
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
# hg38
BiocManager::install("org.Hs.eg.db")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# mm39
BiocManager::install("org.Mm.eg.db")
BiocManager::install("TxDb.Mmusculus.UCSC.mm39.refGene")
```

## Installation

``` r
library(devtools)
devtools::install_github("raulsanzr/methplot")
```

## Usage

### Methylation plot

``` r
meth.plot(genome="hg19", chr="chr11", start=27015473, end=27015991, sites=CpGs, regions=DMR.list, group=metadata$Condition)
```

![](https://github.com/raulsanzr/methylation/blob/main/docs/refs/DMR_1.png)<!-- -->


## Author

Raúl Sanz (2022)
