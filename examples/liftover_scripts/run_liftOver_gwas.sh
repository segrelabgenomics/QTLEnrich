#!/bin/bash

# This Python script runs liftOver_gwas.py to liftOver an inputted GWAS from genome build GRCh37 to GRCh38
# Author: Andrew Hamel
# Date: 18 February 2020
# Written in the Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA
#
# Written in Python 3
# Requires Pandas
# Utilizes subprocess to run liftOver
#

# example using Height GWAS computed from the UK Biobank by the Neale Lab (http://www.nealelab.is/uk-biobank)

gwas="../50_irnt.gwas.imputed_v3.both_sexes.tsv" # file downloaded from: http://www.nealelab.is/uk-biobank, and saved in the examples/ directory.
chain_file="../hg19ToHg38.over.chain.gz" # liftOver chain file from genome build GRCh37 to GRCh38 downloaded from: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/
output="../50_irnt.gwas.imputed_v3.both_sexes.GRCh38.tsv.gz"

#purge and load modules as needed
# load liftOver script downloaded from UCSC: http://hgdownload.soe.ucsc.edu/admin/exe/ and point to appropriate path on local computer/server

module purge
module load <>

./liftOver_gwas.py --gwas_file $gwas \
                   --chain_file $chain_file \
                   --output_file $output

#removes mapped and unmapped files generated during liftOver
rm mapped unmapped


