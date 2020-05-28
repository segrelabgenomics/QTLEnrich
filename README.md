# QTLEnrich v2

**This repository contains the software implementation of QTLEnrich (v2)**, a tool that assesses enrichment of common variant associations for a given complex disease or trait among molecular quantitative trait loci (QTLs), such as expression QTLs (eQTLs) and splicing QTLs (sQTLs) in a given tissue, correcting for potential confounding factors (MAF, distance to TSS, local LD).

**QTLEnrich (v2)** is an extension of eQTLEnrich published in Gamazon*, Segre*, van de Bunt*, et al., Nature Genetics (2018), 50(7):956-967, that accounts for increased discovery power of eQTLs in GTEx release v8 (Genome build hg38) and additional QTL types. The updated version of QTLEnrich is described in the GTEx v8 main paper: [GTEx consortium, BioRxiv 2019](https://www.biorxiv.org/content/10.1101/787903v1) (provisionally accepted at Science 2020), where we applied QTLEnrich to GWAS of 87 complex diseases and traits and eQTLs and sQTLs in 49 tissues. A methods paper is also under preparation.

**Authors**: Andrew Hamel, John Rouhana, written in the Segre lab under the guidance of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA

For questions or comments regarding this tool, please contact Andrew Hamel at andrew\_hamel at meei dot harvard dot edu, and Ayellet Segre at ayellet\_segre at meei dot harvard dot edu.

Date: January 8, 2020

## Repository structure

``src``: the directory contains scripts for the software pipeline and for plotting with relevant instructions

``data``: the directory contains input files required to run QTLEnrich for GTEx v8 QTLs

``examples``: sample data for running QTLEnrich on a given trait and tissue pair

## Dependencies

QTLEnrich was written in Python (at least 3.6) and R (at least 3.6.1) and requires the following libraries and modules:

  * [pandas 0.25.0](https://pandas.pydata.org/)
  * [numpy 1.16.2](https://www.scipy.org/install.html)
  * [scipy 0.17.1](https://www.scipy.org/install.html)
  * [matplotlib 2.2.2](https://matplotlib.org)
  * [ggplot2 3.2.1](https://ggplot2.tidyverse.org)
  * [dplyr 0.8.3](https://dplyr.tidyverse.org)
  * [readr 1.3.1](https://readr.tidyverse.org)
  * [cowplot 1.0.0](https://wilkelab.org/cowplot/index.html)

## Step by step description for running QTLEnrich

Step by step description for running QTLEnrich, preprocessing of input files, and plotting can be found in [src/README.md](https://github.com/segrelabgenomics/QTLEnrich/src/README.md). Instructions are centered on using GTEx v8 data, but can be applied to any non-GTEx QTL datasets.

1. Download and prepare format of QTL results files
2. Compute LD proxy variants for variants tested in QTL study using PLINK
3. Prepare null variants and confounders tables
4. Prepare GWAS input file/s for analysis
5. Run QTLEnrich
6. Generate plots as desired

## License

The QTLEnrich software is distributed under the terms of the BSD 3-Clause License. See [LICENSE.txt](https://github.com/segrelabgenomics/QTLEnrich/LICENSE.txt) file for more details.

## Citation

GTEx Consortium, "The GTEx Consortium atlas of genetic regulatory effects across human tissues", BioRxiv 2019 (accepted at Science 2020) (https://www.biorxiv.org/content/10.1101/787903v1)

Methods preprint in preparation, coming soon!

Last updated: May 26, 2020

