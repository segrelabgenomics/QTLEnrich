# QTLEnrich v2
Assessing enrichment of complex disease or trait associations among QTLs

This repository contains the software implementation of QTLEnrich (v2), which is a tool used to assess enrichment of GWAS variant associations for a given complex disease or trait among molecular quantitative trait loci (QTLs), such as expression QTLs (eQTLs) or splicing QTLs (sQTLs) in a given tissue, correcting for potential confounding factors (MAF, distance to TSS, local LD).

QTLEnrich (v2) is an extension of eQTLEnrich published in Gamazon, Segre et al., Nature Genetics (2018), 50(7):956-967, adapted for the increased discovery power in GTEx release v8 (genome build hg38) and additional QTL types. 

Authors: Andrew Hamel, John Rouhana, written in the Segre lab under the supervision of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA

For questions or comments regarding this tool, please contact Andrew Hamel,andrew\_hamel at meei dot harvard dot edu, and Ayellet Segre at ayellet\_segre at meei dot harvard dot edu.

## Repository structure

`src`: the directory contains scripts for the software pipeline and for plotting with relevant instructions

`data`: the directory contains input files required to run QTLEnrich for GTEx v8 QTLs

`examples`: sample data for running QTLEnrich on a given trait and tissue pair


## Dependencies 

QTLEnrich was written in Python (at least 3.6) and R (at least 3.6.1) and requires the following libraries and modules:

  * [pandas 0.25.0](https://pandas.pydata.org/)
  * [numpy 1.16.2](https://www.scipy.org/install.html)
  * [scipy 0.17.1](https://www.scipy.org/install.html)
  * [matplotlib 2.2.2](https://matplotlib.org)
  * [ggplot2 3.2.1](https://ggplot2.tidyverse.org)
  * [dplyr 0.8.3](https://dplyr.tidyverse.org)
  * [readr 1.3.1](https://readr.tidyverse.org)

## Step by step description of running QTLEnrich

Step by step description of running QTLEnrich, preprocessing input files, and plotting can be found in [src/README.md].

COMING VERY SOON!

## License

The QTLEnrich software is distributed under the terms of the BSD 3-Clause License. See [LICENSE] (license.txt) file for more details.

## Citation

Preprint coming soon on BioRxiv!

Last updated: January 21, 2020
