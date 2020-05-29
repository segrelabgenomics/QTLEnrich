# QTLEnrich v2
 
**This repository contains the scripts and instructions for running QTLEnrich, a tool used to assess enrichment of complex disease or trait associations among expression quantitative trait loci (eQTLs) or splicing QTLs (sQTLs) in a given tissue, correcting for potential confounding factors (MAF, distance to TSS, local LD).**

**QTLEnrich (v2)** is an extension of eQTLEnrich published in Gamazon*, Segre* *et al., Nature Genetics (2018), 50(7):956-967, modified to account for the increased discovery power of eQTLs in GTEx release v8 (genome build hg38) and other QTL types. Updates of QTLEnrich v2 are described in the GTEx v8 main paper: [GTEx consortium, BioRxiv 2019](https://www.biorxiv.org/content/10.1101/787903v1) (accepted to Science 2020), where we applied **QTLEnrich** to GWAS of 87 complex diseases and traits and eQTLs and sQTLs in 49 tissues. A methods paper is also under preparation.

**QTLEnrich** can also output a list of target genes of eQTLs or sQTLs with top ranked GWAS p-values (e.g., p<0.05) that can be inputted into **GeneEnrich** to test for over-representation of trait-associated target genes in specific biological pathways or processes, if enrichment of GWAS associations among QTLs is found. **GeneEnrich** software will be added to Github soon.

**We provide relevant confounder and null files for GTEx v8 eQTLs, sQTLs, and independent eQTLs**, that can be downloaded from the GTEx portal: https://www.gtexportal.org/home/datasets or from https://github.com/broadinstitute/gtex-v8. 

**We also provide scripts to process QTLEnrich input files for non-GTEx QTL and any GWAS study.**

**Authors:** Andrew Hamel, John Rouhana, written in the lab of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA

**Date:** October 2, 2019

For questions or comments regarding this tool, please contact Andrew Hamel: andrew\_hamel at meei dot harvard dot edu and Ayellet Segre at ayellet\_segre at meei dot harvard dot edu.

## Dependencies 

QTLEnrich was written in Python (at least 3.6) and requires the following libraries and modules

  * [pandas 0.24.2](https://pandas.pydata.org/)
  * [numpy 1.16.2](https://www.scipy.org/install.html)
  * [scipy 0.17.1](https://www.scipy.org/install.html)
  * See additional libraries required for plotting options below.

## Running Environment and Run Time

We recommend to reserve ~75GB of memory when running **QTLEnrich**.

These are estimated run times for a sample GWAS (UK Biobank Height GWAS) and eQTLs from GTEx tissues using a default of up to 100,000 permutations:

 * ~1,000 significant eQTLs: Kidney Cortex: 12-15 minutes
 * ~10,000 significant eQTLs: Heart Atrial Appendage: 15-20 minutes
 * eQTL sets from 49 GTEx v8 tissues: 300 minutes (5 hours)

## Repository structure

``src``: the directory contains scripts for the software pipeline and for plotting with relevant instructions

``data``: the directory contains input files required to run QTLEnrich for GTEx v8 eQTLs and sQTLs

``examples``: sample data for running QTLEnrich on a given trait and tissue pair

# Instructions on how to run QTLEnrich.py

### To run QTLEnrich you have to complete the following steps described in detail below:

A sample command for running QTLEnrich.py on Standing Height GWAS from the UK Biobank and Whole blood QTLs from GTEx v8 is given in this shell script:
**../examples/example\_qtlenrich\_command.sh**

All required precomputed input files for running **QTLEnrich** on GTEx v8 QTLs, sQTLs, and independent eQTLs (as described in steps 2-6) are provided in the *data* directory.

### A. Data preparation and preprocessing

1. Download QTL results files and check in proper format
2. Prepare Null variants table per tissue
3. Compute LD proxy variants for all variants tested in QTL study
4. Process GENCODE file
5. Generate Confounders Table
6. Download and preprocess GWAS summary statistics results file

### B. Run *QTLEnrich*

Once all input files are ready, the QTL-GWAS enrichment analysis can be run using **QTLEnrichV2.py**. A sample command with input and output variables is provided here: **sample_qtlenrich_command.sh**. See below for a detailed description of all input and output variables and files.

### C. Plotting options


## A. Data Preparation and Pre-Processing

The input files available for **QTLEnrich** are based on eQTL and sQTL results from 49 tissues in GTEx release v8. To run **QTLEnrich** on other QTL studies, you will need to perform steps 1-4 described below, before running QTLEnrichV2.py.

### 1. Download and store eQTL or sQTL results files in the *data/* directory in proper format.

**To run QTLEnrich on GTEx v8**, you will need to download the files with significant eQTL and/or sQTLs from the GTEx portal: https://www.gtexportal.org/home/datasets or from GitHub: https://github.com/broadinstitute/gtex-v8, and save them in the *data*directory in your clone of this repository. *QTLEnrich* can be run on all of the following QTL types and files with no processing needed (if inputted in their original format):

* eQTLs, including the most significant eQTL variant per eGene: GTEx_Analysis_v8_eQTL.txt.gz retrieved from GTEx_Analysis_v8_eQTL.tar
* sQTLs, including most significant sQTL variant per sGene: GTEx_Analysis_v8_sQTL.txt.gz retrieved from GTEx_Analysis_v8_sQTL.tar
* Independent eQTL signals: GTEx_Analysis_v8_eQTL_independent.tar
* Independent sQTL signals: GTEx_Analysis_v8_sQTL_independent.tar
* Fine-mapped eQTLs and sQTLs: GTEx_v8_finemapping_DAPG.tar
* Cell type interaction eQTLs: GTEx_Analysis_v8_ieQTL.tar
* Cell type interaction sQTLs: GTEx_Analysis_v8_isQTL.tar
* Trans eQTLs: GTEx_Analysis_v8_trans_eGenes_fdr05.txt
* Trans sQTLs: GTEx_Analysis_v8_trans_sGenes_fdr05.txt

**To run QTLEnrich on your own QTL study**, you will need to prepare a best QTL per gene file (one row per gene) with the following headers: Required headers: variant_id, gene_id. Optional headers: qval (QTL q-value; necessary if list of best QTL per gene is not filtered on q-value cutoff), TSS_Tissue (distance of QTL variant to the transcription start site (TSS) of the target gene). TSS_Tissue can also be computed within QTLEnrich.
The prefix of the QTL file name must contain the tissue name.

All other required input files for analyzing GTEx v8 QTLs, defined below (steps 2-6) are provided in the *data* directory.

### 2. Curation of a Null Table

A null table needs to be generated that defines which variants are non-significant eQTLs or sQTLs. We defined two types of null variant sets - tissue-dependent and global, if multiple tissues are considered. The tissue-dependent null sets contain one column per tissue (listed in header), with 1 representing a null variant, and 0 representing either a significant variant in the given tissue or a variant not tested in that tissue which may be a significant or non-significant QTL in any of the other tissues. The tissue names for the headers are extracted from the prefix of the QTL file names, thus the file name prefixes must contain the tissue name (see GTEx eQTL file names for example). The global null set is defined as all variants tested that are a non-significant variant-gene pair eQTL or sQTL in any of the tissues tested in the study (unique union of all tissue-dependent null variants).

**For GTEx v8 eQTLs and sQTLs**, we defined both tissue-dependent and global null sets for the 49 tissues in GTEx v8: (i) (Default on **QTLEnrich**) Tissue-dependent null sets were defined per tissue as the subset of all variants within +/-1Mb around the TSS of all genes expressed in the given tissue, that are non-significant variant-gene pair eQTLs or sQTLs (FDR>5%) in any of the 49 tissues tested in GTEx v8. We defined tissue-dependent null sets for each of the 49 tissues. (i) The global null set was defined as all tested variants (+/-1Mb around all genes’ TSS) that are a non-significant variant-gene pair eQTL or sQTL (FDR>5%) in all of the 49 tissues tested.

Null table for GTEx v8 **GTEx_v8_Null_Table_49Tiss_Global.txt.gz** can be downloaded from [here](http://segrelab.meei.harvard.edu/software), and should be saved in *data* directory.

**The following scripts and commands were used to generate the null table:**

#### 2.1. Generate file with unique list of variant IDs of all variants tested in the QTL analysis:

```
./extract_unique_variants.py --directory QTL_directory --file_name qtl_file_suffix --output_file unique_list_variants
```

##### Inputs to extract_unique_variants.py to generate unique list of variants tested

* ``--directory``: directory containing qtl files for all tissues tested

* ``--file_name``: the suffix of all QTL tissue files with all the variant-gene pairs tested saved in the directory (e.g., '.allpairs.txt.gz' in GTEx release v8)

* ``--output_file``: name of outputted file with unique list of variants tested, e.g. unique_list_variants (inputted into create_null_variants.py)

We also provide a shell script to run *extract_unique_variants.py*.

```
./run_extract_unique_variants.sh
```

#### 2.2. Generate list of significant QTLs (e.g., eQTLs and sQTLs if both computed)


```
./extract_unique_variants.py --directory QTL_directory --file_name qtl_file_suffix --output_file unique_list_significant_variants_file
```

##### Input variables the same as in 2.1, except for the \-\-file_name:

* ``--directory``: directory containing qtl files for all tissues tested

* ``--file_name``: the suffix of all QTL tissue files with just the significant variant-gene pairs (e.g., '.v8.signif_variant_gene_pairs.txt.gz' for significant eQTL files and '.v8.sqtl_signifpairs.txt.gz' for significant sQTL files

* ``--output_file``: name of outputted file with unique list of significant variants tested, e.g. significant_variants (inputted into create_null_variants.py)

**Note:** if running this on both eQTLs and sQTLs, this command needs to be run on each QTL type independently, and the user needs to generate one concatenated unique list of significant variants.

#### 2.3. Generate list of null variants:

```
./create_null_variants.py --all_variants_file unique_list_variants --significant_variants_file significant_variants --output_file null_variants
```

##### Inputs to generate null variants

* ``--all_variants_file``: File containing unique list of all variants tested (generated in step 2.1).

* ``--significant_variants_file``: File containing list of significant variants (generated in step 2.2).

* ``--output_file``: name of outputted null variant file, e.g. null_variants

We also provide a shell script for to run *create_null_variants.py*.

```
./run_create_null_variants.sh
```

#### 2.4 Generate the null variants table:

```
./Generate_Null_Table.py --qtl_directory directory_for_all_pairs_files --File_Extension .allpairs.txt.gz --Null_Variants_List null_variants --QTL_Significant_Variants_List significant_variants_qtl
```

or run the following shell script:

```
./run_Generate_Null_Table.sh
```

##### Required inputs to generate null table:

* ``--qtl_directory``: directory containing the QTL files with all variant-gene pairs tested (in GTEx v8: *allpairs.txt.gz). Only required column/header: 'variant_id'. 

* ``--File_Extension``: suffix of the QTL files with all variant-gene pairs tested (in GTEx v8: *allpairs.txt.gz) in the directory.

* ``--Null_Variants_List``: unique list of all null variants across all tissues (generated in step 2.3). This file should consist of one column with header "variant\_id". One variant per line.

##### Optional inputs: 

* ``--QTL_Significant_Variants_List``: additional list of significant QTL variants to remove from the null variant set. This option can be useful if an additional QTL type or QTL tissue is computed at a later stage. This file is in the same format as \-\-Null\_Variants\_List.

Output file is a 1/0 tables with one row per variant that specifies if the variant is null (=1) per tissue or globally. This file is inputted into QTLEnrich (see details below).

### 3. Compute LD proxy variants using PLINLK 1.90

#### 3.1. The following command using Plink 1.9 can be used compute number of LD proxy variants.  For the analysis of GTEx v8 QTLs, we used the European subset of the WGS VCF from GTEx release v8 

```
plink --bfile binary_file_name --keep eur_ancestry --r2 --ld-snp-list gtex_variants --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0.5 --out "$output_dir"
```

**Note:** This command was run on variants for each chromosome.

#### Inputs taken from Plink 1.9 documentation: [PLINK documentation] (https://www.cog-genomics.org/plink/1.9)

* ``--bfile``: Prefix of Plink binary fileset. Can run by chromosome.

* ``--keep``: required if subsetting on a particular ethnicity (e.g. European ancestry)

* ``--ld-snp-list``: list of variants for which to find LD proxy variants. Variants must be in the same format as the .bim file in PLINK. The variant format is different by genome build and resource.

* ``--ld-window``: maximum number of variants in LD window

* ``--ld-window-kb``: window size of LD block

* ``--ld-window-r2``: r2 cutoff for LD proxy variants

* ``--out``: output directory


#### 3.2. Parse Plink LD proxy output file into a table that contains number of LD proxy variants per variant that is inputted into *Generate_Confounders_Table.py*:

```
./create_ld_table.py --plink_output_directory output_ld_directory --output_file number_ld_proxy_per_variant_file 
```

#### Inputs for create_ld_table.py

* ``--plink_output_directory``: directory containing outputs from plink command used to compute number of ld proxy variants per variant. If the Plink --r2 command above was run on each chromosome separately, the script *create_ld_table.py* will loop over all output files found in this directory

* ``--output_file``: output file with number of LD proxy variants per variant, e.g. number_ld_proxy_per_variant_file

We also provide a shell script for to run *create_ld_table.py*

```
./run_create_ld_table.sh
```

### 4. Processing the GENCODE file

This step extracts all gene start positions needed to compute the distance to transcription start site (TSS) of all variants tested for eQTL or sQTL analysis.

A GENCODE GTF file in the relevant version needs to be downloaded from https://www.gencodegenes.org/human/.
The GENCODE GTF file is processed into a table-like document to make it more readably in pandas in the script that generates the Confounders table: *Generate_Confounders_Table.py*, using the gtf parser provided by Kamil Slowikowksi: [https://gist.github.com/slowkow/8101481] (https://gist.github.com/slowkow/8101481)

**For GTEx release v8 QTLs,** we reformatted the GENCODE v26 file download from: https://www.gencodegenes.org/human/release_26.html, using the following command:

```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz
gunzip gencode.v26.annotation.gtf.gz
```

The processed GENCODE v26 file (GENCODE_26_df.txt) is saved in the *data* folder.

### 5. Generating the Confounders Table

The confounders table includes annotations of 3 confounding factors for all variants (null and significant) that were tested for eQTL and sQTL significance (+/-1MB around TSS of all expressed genes), computed for each of the 49 tissues in GTEx v8 independently and globally. Confounding factors include minor allele frequency (MAF), minimum distance to TSS of all genes expressed in the given tissue (tissue-dependent) or across all tissues (global), and number of LD proxy variants at r2>0.5 (LD\_proxy). For the significant eQTLs or sQTLs, distance to TSS is re-computed in **QTLEnrich.py** considering only the eGenes with a significant QTL in the given tissue. MAF and distance to TSS were computed for each tissue, as they vary by sample set and by the set of genes expressed in each tissue (tissue name added to column header). TSS is computed as the variant position minus the TSS position in base pair units, and strand orientation of each gene was accounted for by multiplying the difference by -1 if the gene is on the negative strand (variant upstream of TSS: negative value; variant downstream of TSS: positive value).

**We generated this confounders table for GTEx release v8 eQTLs and sQTLs** using GENCODE v26, saved in *data*. Number of LD proxy variants was computed using r2>0.5 and the European subset of samples in the WGS VCF of GTEx v8, and MAF was taken from the WGS VCF computed considering all 838 donors.

Confounders table for GTEx v8 **GTEx_v8_Confounders_Table_49Tiss_Global.txt.gz** can be downloaded from [here](http://segrelab.meei.harvard.edu/software) and should be saved in *data* directory.

**The following script and command can be used to generate the confounders table for a new QTL dataset.** Pre-requisites includes computing number of LD proxy variants in Step 3, processing the GENCODE GTF file in step 4, and inputting a table with all MAF for each variant. Also, the tissue names will be taken from the prefix of the QTL file names, hence the prefix of the QTL file name must contain the tissue name.

The list of full variants was computed using the script *extract_unique_variants.py*. See above in curation of null table for more information.

``` 
./Generate_Confounders_Table.py --qtl_directory directory_for_all_pairs_files --File_Extension .allpairs.txt.gz --Variants_List unique_list_variants --LD_Proxy 
number_ld_proxy_per_variant_file --GENCODE_file gencode_processed_file --Output_File output_confounder_file 
```

We also provide a shell script to run **Generate_Confounders_Table.py**:
```
./run_Generate_Confounders_Table.sh
```

#### Inputs into Generate_Confounders_Table.py:

* ``--qtl_directory``: directory containing the eQTL files with all variant-gene pairs tested (in GTEx: *allpairs.txt.gz). Three required columns/headers: 'variant_id', 'gene_id', and 'maf'.

* ``--File_Extension``: suffix of the eQTL files of all variant-gene pairs tested (in GTEx: *allpairs.txt.gz) in the directory.

* ``--Variants_List``: Unique list of all variants tested in the QTL analysis (significant and null variants, generated in step 2.1). This file should consist of one column with header "variant\_id". One variant per line.

* ``--LD_Proxy``: table of all variants and the number of LD proxy variants per variant generated in step 3.2. This file should be a tab-delimited table with two columns: 'variant' and 'LD\_Proxy' obtained from *create_ld_table.py*

* ``--GENCODE_file``: gencode dataframe generated in step 4 (provided for GENCODE v26 in this QTLEnrich package). Note: This confounders table was constructed using ENSEMBL gene IDs.

* ``--Output_File``: name of output confounder table. This script will compress the file to save memory. It is therefore recommended to name your file with a .gz extension, e.g. confounders\_table.txt.gz

### 6. Processing the GWAS input file

The GWAS input file can be tab- or space-delimited and requires three headers: "chr", "pos", "gwas\_p\_value", with two optional headers: "effect\_allele" and "non\_effect\_allele". The "effect\_allele" in the GWAS study should match the "effect\_allele" in the QTL study. If using the GTEx variant ID nomenclature, the non\_effect\_allele should match the ref allele, while the effect\_allele should match the alt allele. The file can be saved as .txt or txt.gz.

    The genome build of the GWAS has to match that of the QTL study. For genome build 38, the chr numbers will be listed as "chr" plus the chromosome number, e.g. chr10. If the GWAS file contains variants from chromosomes X,Y, and the mitochondrion (MT), variants must be listed as chrX, chrY, and chrMT, respectively.

        chr    pos	non_effect_allele  effect_allele  gwas_p_value
        chr1   14981    A    G    0.003
        chrX   18200    C    G    0.42
        chr22  50801268 C    T    0.05

    Additionally, it is recommended that the GWAS file name contain the GWAS trait name as it will be used in the output file names, e.g.:

        type_2_diabetes.txt.gz
        coronary_artery_disease.txt.gz

## B. Command for running QTLEnrich:

Command for running **QTLEnrich** with appropriate input and output variables and files is provided in the shell script: **sample_qtlenrich_command.sh**.

## Inputs into QTLEnrichV2.py script:

### REQUIRED INPUTS:

* ``--gwas_file``: The GWAS input file can be tab or space-delimited and requires three headers: "chr", "pos", "gwas\_p\_value", with two optional headers: "effect\_allele" and "non\_effect\_allele". See more details above in Point 6.

    For genomic build GRCh38, the chr numbers will be listed as "chr" plus the chromosome number, e.g. chr10. If the GWAS file contains variants from chromosomes X,Y, and the mitochondrion (MT), variants must be listed as chrX, chrY, and chrMT, respectively. For genome build GRCh37, the chr numbers will be without the prefix 'chr'.

        chr    pos      non_effect_allele  effect_allele  gwas_p_value
        chr1   14981    A    G    0.003
        chrX   18200    C    G    0.42
        chr22  50801268 C    T    0.05

    Additionally, it is recommended the GWAS file name contain the GWAS trait name as it will be used in the output file names, for example:

        type_2_diabetes.txt.gz
        coronary_artery_disease.txt.gz

* ``--qtl_directory``: Directory containing significant QTL results files (e.g., files with most significant QTL variant per gene). Files need to be tab-delimited. QTL files must consist of a column with a variant header "variant\_id" that contains variant IDs in concordance with GTEx nomenclature and "gene\_id" consisting of the ENSEMBL ID. Format of variant IDs should follow the following form: chrN\_pos\_ref\_alt\_GenomeBuild.

Depending on the QTL type, additional columns may be necessary.  When assessing enrichment of eQTLs or sQTLs, a q-value header, "qval" can be added.  Significant QTLs will be filtered on a desired q-value cutoff (See --qtl\_q\_value below for more information).  If this header is not provided, it is assumed that all variants in the file are significant.  Similarly, if assessing enrichment on independent QTLs, an independent ranking header can be added.  These values will be treated as integers with the header of "rank".

        variant_id    gene_id       
        chr15_90160604_G_A_b38    ENSG0000123.4    
        chr3_51671175_A_G_b38     ENSG00003244.17
        chr1_153541765_A_C_b38    ENSG000095.1

A sample input file is provided here: */data/Whole\_Blood.v8.egenes.txt.gz*.

* ``--file_name``: should contain only the suffix of all QTL files that will be read in the QTL Directory, used on the egenes files.

    **List of suffix options:

    * .v8.egenes.txt.gz
    * .v8.sqtl_signifpairs.txt.gz
    * .v8.independent_eqtls.txt.gz

    **Note:** The prefix of each file read with this suffix in the QTL\_Directory must contain the tissue name. For analysis of GTEx QTLs, please use the GTEx nomenclature listed in "Input\_files/GTEx\_Tissue\_Names.txt". For example, in the file name: "Whole\_Blood.variants.txt.gz", the tissue "Whole\_Blood" is the prefix, and ".variants.txt.gz" is the suffix. The tissue name will be used for the "tissue\_dependent" option in the null variant selection if chosen (see below)

* ``--qtl_type``: Denotes QTL type analyzed in **QTLEnrich**.  Current options include best\_eqtl, best\_sqtl, conditional\_independent, dapg\_independent, ieqtl, isqtl.

* ``--confounders_table``: The confounders table will include a list of all unique variants (null and significant) and their respective confounding factors for each of the tissues included in the QTL study, e.g., 49 tissues in GTEx release v8 ("tissue_dependent"), and considering genes expressed in at least one of the tissues ("global"). Confounding factors include minor allele frequency (MAF), distance to transcription start site (TSS), and number of LD proxy variants at r2>0.5 (LD\_proxy). MAF and distance to TSS are computed per tissue (tissue and name added to column header) and across all tissues ("Global", in which case minimum distance to TSS across all tissues is taken). Strand orientation was taken into consideration in calculating the distance to TSS, given in base pair units.

We precomputed this file for GTEx release v8: **GTEx_v8_Confounders_Table_49Tiss_Global.txt.gz**, which can be downloaded from [here](http://segrelab.meei.harvard.edu/software) and should be saved in *data* directory.

        variant_id              MAF_Whole_Blood TSS_Whole_Blood MAF_Liver  TSS_Liver      LD_proxy
        chr15_90160604_G_A_b38  0.272556        -57             0.282      123        7
        chr3_51671175_A_G_b38   0.147556        5               0.342     -90         0
        chr1_153541765_A_C_b38  0.101504        100             0.123      87         5
          

* ``--null_table``: Table that contains all variants tested in the QTL anlaysis and specifies which are null per tissue or globally. For example, the first variant below is a null variant only in Whole Blood, the second variant is a null variant only in Liver, and the third variant is null in both tissues. See Curation of Null Table in Section 2 above for more information.  We precomputed this file for GTEx release v8: GTEx_v8_Null_Table_49Tiss_Global.txt.gz, which can be downloaded from [here](http://segrelab.meei.harvard.edu/software), and should be saved in *data* directory.

        variant_id         Tissue_Whole_Blood Tissue_Liver  Global   
        chr15_90160604_G_A_b38  1        0               1          
        chr3_51671175_A_G_b38   0        1               1
        chr1_153541765_A_C_b38  1        1               1

### OPTIONAL INPUTS:

* ``--subset_genes``: indicate if user wishes to subset on protein coding and lincRNA genes. Requires \-\-gencode\_file (default: True).

* ``--keep_gene_id_suffix``: indicate if user wishes to keep decimal version suffix on ENSEMBL gene ids (default: True).

* ``--compute_tss_distance``: indicate if user wishes QTLEnrich to compute TSS distance for significant QTLs. Requires \-\-gencode\_file (default: True).

* ``--gencode_file``: gencode_file used for computing TSS distance for significant QTLs and their target gene, and for subsetting on protein coding and lincRNAs genes. This file is provided for GENCODE v26 which is compatible with GTEx v8 (/data/GENCODE_26_df.txt).

* ``--null_option``: indicates if null sampling is performed using a tissue-dependent null or global null. Options include "tissue\_dependent" or "global" (default: "tissue\_dependent").

* ``--genome_build``: human genome build of the variant chromosome positions used.  If analyzing eQTLs/sQTLs from GTEx release v8, variants must be in genome build GRCh38 (hg38); for GTEx release v7, variants are in GRCh37 (hg19).  Importantly, the GWAS variants and QTLs must be from the same genome build. Input format: "b38" or "b37". Default is "b38".

* --qtl_q_value: q-value cutoff if columnn with eGene or sGene q-value is provided. Default cutoff is 0.05.

* ``--gwas_p_value``: GWAS p-value cutoff used for QTL-GWAS enrichment analysis. Default cutoff is 0.05.

* ``--independent_ranking``: denotes rank of independent QTLs to be used in enrichment analysis (e.g., 1 for primary signal, 2 for secondary signal, etc.), if chosen QTL type is 'independent'. Default ranking is 1.

* ``--lambda_factor``: denotes cutoff for lambda value when computing our empirical Pi1. Default lambda value is 0.8.

* ``--lower_bound_permutations``: designates the minimum number of permutations the user would like to perform. Default is 1,000 permutations.

* ``--upper_bound_permutations``: designates the maximum number of permutations that will be performed under our adpative permutation scheme. Default is 100,000 permutations.

* ``--num_quantiles``: number of quantile bins used for each confounding factor when matching on confounders with the significant QTLs in the permutation scheme. Default is 10.

* ``--tissue_name``: this file contains a list of tissues names to be analyzed, in case the user wants to analyze the significant QTLs of only a subset of tissues. Each row should contain one tissue, and tissue names should be identical to the tissue names used in the Confounders table and Null table. If analyzing GTEx QTLs, tissue names are listed in "../data/GTEx\_Tissue\_Names.txt". For example:

        Whole_Blood
        Nerve_Tibial
        Adrenal_Gland
        Brain_Spinal_cord_cervical_c-1

*  ``--output_directory``: option to provide directory path or create a new output directory where all output files will be saved. Default is to create a new output directory (e.g., Output\_QTLEnrich\_exp\_label\_DATE).

*  ``--exp_label``: unique name or details on the analysis run, which will be added to the output directory and output file names.

### Options for generating GeneEnrich input files:

* ``--GeneEnrich_Input``: option to generate input files for GeneEnrich (no value needed).

These two options are needed only if the 'independent QTL' type is chosen:

* ``--eGenes_directory``: directory that contains the eGene or sGene results files (list of all genes tested with the most significant QTL per gene and gene level q-value)
This is need to extract the list of non-significant genes for the independent QTLs.

* ``--eGenes_file_name``: File name suffix of the eGenes or sGenes results files in the \-\-eGenes\_directory. Similar to \-\-file\_name option for main set of QTLs.
Required columns: "gene\_id", "variant\_id", "qval"

## Output Tables Format

**QTLEnrich produces the following output files.** A description of the column headers of all files can be found in: **../qtlenrich_output_header.txt**

1. Main results table with the QTL-GWAS enrichment p-values and adjusted fold-enrichment (*QTLEnrich_{qtl_type}_{trait_name}_{exp_label}_{date}.tsv*)
2. GWAS p-values of the set of significant QTLs (*QTLEnrich_significant_gwas_p_values_{qtl_type}_{trait_name}_{tissue}_{date}.tsv*)
3. Variant IDs of the first N (default 1000) permuted null variants matched on three potential confounding factors relative to the significant QTL variants (*QTLEnrich_null_variant_id_matrix_{qtl_type}_{trait_name}_{tissue}_{date}.tsv*). This table can be used for replication analysis in a second independent GWAS.
4. GWAS p-values of N (default 1000) matched null variants in table above (in point 3) (*QTLEnrich_matched_null_gwas_p_values_{qtl_type}_{trait_name}_{tissue}_{date}.tsv*).
5. _GeneEnrich_ Input Files: if the \-\GeneEnrich\_Input option is indicated, two GeneEnrich input files are generated: (i) table of significant eGenes or sGenes (e.g., q<0.05) with a best QTL variant or an independent QTL at a 
specified QTL rank (default 1) (depending on specified QTL type) with GWAS p-value below a given cutoff (default: GWAS P<0.05). Columns contain: 'ENSEMBL ID', 'variant_id' (of best QTL per gene); (ii) null set of genes: list of 
all genes (genes with no significnat and non-significant eGenes or sGenes) with QTL whose GWAS p-value is above a given cutoff (default: GWAS P>0.05). One column: 'ENSEMBL ID'. If 'best QTL' option is chosen, the GWAS p-value of 
the best QTL variants for all genes is evaluated. If 'independent QTL' option is chosenn the null set of genes is defined as the significant eGenes or sGenes whose independent QTL at a pre-specified rank has a GWAS p-value above a 
given cutoff plus the non-significant genes (e.g., q>0.05) whose best QTL has a GWAS p-value above a given cutoff. Examples of these output files can be found in the *../examples/Output_QTLEnrich_example_Mar02_2020/* folder.

# C. Plotting Functions of QTLEnrich Output Data

**This repository contains the scripts needed to visualize results and outputs from *QTLEnrich* (v2). The possible plots include:**

1. Q-Q plots to inspect the GWAS p-value distributions for a given set of tissue-trait pairs or QTL types.
2. Forest plots of the adjusted fold-enrichment for multiple tissue-trait pairs and QTL types.
3. Bar plot of the estimated number of trait associations among QTLs based on adjusted fold-enrichment or empirical Pi1.
4. Box or violin plots of the adjusted fold-enrichment, number of trait associations, or empirical Pi1 for given trait and multiple tissues and/or QTL types.

**The QQ-Plot module was written in Python 3.6+ requires the following libraries:**

  * [pandas 0.25.0](https://pandas.pydata.org/)
  * [numpy 1.16.2](https://www.scipy.org/install.html)
  * [scipy 0.17.1](https://www.scipy.org/install.html)
  * [matplotlib 2.2.2](https://matplotlib.org)

**The forest plots, bar plots, and box plots modules were written in R 3.6.1 and requires the following packages:**

  * [ggplot2 3.2.1](https://ggplot2.tidyverse.org)
  * [dplyr 0.8.3](https://dplyr.tidyverse.org)
  * [readr 1.3.1](https://readr.tidyverse.org)
  * [cowplot 1.0.0](https://wilkelab.org/cowplot/index.html)

### Sample commands for different plotting options

    ./sample_qqplot_command.sh

    ./sample_forestplot_barplot_command.sh

    ./sample_boxplot_violinplot_command.sh

We generated sample figures for Height GWAS and Whole Blood eQTLs and sQTLs in: **../examples/Output_QTLEnrich_example_Mar02_2020/**

## 1. Q-Q Plots

Example generation of Q-Q plots on Height GWAS and Whole Blood eQTLs is given in: **./sample_qqplot_command.sh** that is a wrapper of **./qqplot_qtlenrich.py**.

### Required Inputs

* ``--qtl_file (outputted from QTLEnrich)``:  GWAS p-values of set of significant QTLs. Tab-delimited file comprised of two columns: 'variant' and 'gwas_p_value'. For example:

        variant                   gwas_p_value
        chr15_90160604_G_A_b38    0.003
        chr1_153541765_A_C_b38    0.05

### Optional inputs

**For the background**, user can choose to plot the GWAS p-values of the full GWAS or of null variants matched on confounding factors compared to significant QTLs (see next two options):   

1. ``--gwas_background``: Option to plot a Q-Q plot with all GWAS variant p-values as background. This plot consists of the full GWAS and significant QTLs.

If this option is chosen, define GWAS input file here:

* ``--gwas_file``: GWAS p-values of all genome-wide variants in given GWAS; input file must include a 'gwas_p_value' column (can be a single column).

        gwas_p_value
        0.0034
        0.513
        0.00001

2. ``--matched_null_background``: Option to plot Q-Q plot with a matched-null background. This plot overlays GWAS p-values for N sets of matched null variants with the GWAS p-values of the significant QTL set.

If this option is chosen, define null variant input file here:

* ``--pvalue_null_matrix``: Matrix outputted from **QTLEnrich** with GWAS p-values of confounder-matched null variants for multiple (default=1000) permutations (rows=variants, columns: permutations).

    Example file: "../examples/Output_QTLEnrich_example_Mar02_2020/output_files/QTLEnrich_matched_null_gwas_p_values_best_eqtl_50_irnt.gwas.imputed_v3.both_sexes.GRCh38_Whole_Blood_Mar02_2020.txt"

### Other options:

* ``--num_permutations``: Indicates number of matched null permutations whose GWAS p-values will be plotted as background on the Q-Q plot (default: 100).

* ``--confidence_interval``: Option to shade confidence interval of GWAS p-values of matched null variants instead of plotting GWAS p-values of N sets of matched null variants across \-\-num\_permutations (default: False).
If option is chosen, user can define the lower and upper bound confidence intervals (default: 95% CI) with the two options below:

    * ``--upper_bound_ci``: Select upper bound confidence interval (default: 97.5 percentile)

    * ``--lower_bound_ci``: Select lower bound confidence interval (default: 2.5 percentile)

* ``--upper_bound_clip``: Set upper bound p-value on Q-Q plot (default: no clip, type is an integer,e.g. 20)

* ``--trait_label``: Label for GWAS trait (default : Trait)

* ``--qtl_type``: Type of QTLs in the plot, e.g., 'eQTL', 'sQTL' (default: Significant QTLs)

* ``--tissue_label``: Name of tissue displayed in the plot (default: Tissue)

* ``--show_title``: Option to include title in plot (default: False)

    * ``--qqplot_title``: Option to input own title, requires \-\-show\_title set to True (default: QQ-Plot)

* ``--output_directory``: Name of output directory for plots (default: Same directory name as that specified for the QTLEnrich results files, e.g., Output\_QTLEnrich\_exp\_label\_DATE).

* ``--exp\_label``: Add label to filenames (default: no label)

* ``--x_inches``: Number of horizontal inches on plot (default: 12)

* ``--y_inches``: Number of vertical inches (default: 10)

* ``--dpi``: Display pixels per inch (default: 100)


## 2-3. Forest Plots and Bar Plots

Example for generating a forest plot of the adjusted fold-enrichment of Height GWAS associations among Whole blood eQTLs and sQTLs side by side with a bar plot of estimated number of height associations amongt blood eQTLs based on Pi1 is given in **./sample_forestplot_barplot_command.sh** that is a wrapper of **./forest_boxplot_module.R**.

#### Required Inputs

* ``--qtlenrich_file1``: QTLEnrich results output file used as input into all types of plotting features (QTLEnrich_{qtl_type}_{trait_name}_{exp_label}_{date}.tsv). 

**This ./forest_boxplot_module.R module requires at least one of the following plotting options:**

* ``--forest_plot``: Option to generate forest plot of the adjustment fold-enrichment (mean and 95% confidence interval) per tissue.

* ``--bar_plot``: Generates bar plot. If bar\plot option is chosen, user can use the following option to choose the summary statistic to display:
		
	* ``--bar_box_violin_option``: Options include: 'pi1\_trait', which plots number of trait associations based on an empirical Pi1; 'pi1', which plots the empirical Pi1 rate; 'adjusted\_fold', which plots adjusted fold-enrichment (default: pi1\_trait)

### Optional Inputs

If user is comparing two QTL outputs from **QTLEnrich** (e.g., eQTL and sQTL):

* ``--qtlenrich_file2``: Second QTLEnrich output file if comparing results from two QTL outputs (e.g., eQTL and sQTL) (QTLEnrich_{qtl_type}_{trait_name}_{exp_label}_{date}.tsv).

* ``--paired_plot``: Indicates user wants to generate a paired forest plot or bar plot (e.g., two QTL types).

### Other options:

* ``--label_file1``: Label for first QTLEnrich output file (e.g., eQTL). If label is more than one word, e.g. Type 2 Diabetes, please add quotes to variable or string in shell script (default: "QTL type 1"):

```
./forest_boxplot_module.R --label_file1 "Type 2 Diabetes"
```

* ``--label_file2``: Label for second QTLEnrich output file (e.g., sQTL). For format of input, please see \-\-label\_file1 (default: "QTL type 2") 

* ``--significant_option``: Flag specifies generating a plot including only the significant tissues. The significance cutoff is defined in following option:

	* ``--identify_significant_cutoff``: Enrichment p-value cutoff for significant tissues (e.g. Enrichemnt P<0.001). If not indicated, default is Bonferroni cutoff (default: 0.05/number of tissues in QTLEnrich output file).

* ``--identify_number_significant_tissues``: If user wants to plot a subset of top ranked tissues, ranked based on their adjusted fold-enrichment, e.g. first 10 tissues.

* ``--remove_subset_tissues``: To exclude specific tissues, include text file with list of tissues to exclude. This file is comprised of one column with header 'Tissue (QTL)' and tissue names should match input QTL files.
Official GTEx tissue names are listed in: "../data/GTEx\_Tissue\_Names.txt".

        Tissue (QTL)
        Whole_Blood
        Nerve_Tibial
        Brain_Spinal_cord_cervical_c-1

* ``--gtex_tissue_names``: If QTLs are from GTEx, this option uses GTEx tissue names as tissue labels. Tissue names can be found in this file downloaded from the GTEx portal for v8: "../data/gtex_tissue_colors_v8.txt."

* ``--side_by_side``: Option to plot forest plot and bar plot panels vertically side by side (uses cowplot). If option is indicated, it requires both \-\-forest\_plot and \-\-bar\_plot options to be specified. If both the \-\-forest\_plot and \-\-bar\_box\_violin\_option options are specified, without the \-\-side\_by\_side option, two separate plots are generated one for the forest plot of the adjusted fold-enrichment and the other for the bar plot of the summary statistic specified here. 

* ``--output_directory``: Place plots in directory of choosing (default: Same directory name as that specified for the QTLEnrich results files, e.g., Output\_QTLEnrich\_exp\_label\_DATE)

* ``--exp_label``: Add label to filenames (default: no label).

## 4. Box Plots or Violin Plots

This plotting feature is relevant when you test QTLs in multiple tissues per trait or vice versa.

Example of plotting side by side box plots of adjusted fold-enrichment of Height associations among GTEx eQTLs and sQTLs from 49 tissues is given in **./sample_boxplot_violinplot_command.sh** which is a wrapper of **./forest_boxplot_module.R**. Similar plots can be generated for estimated number of trait associations or empirical Pi1 amonng QTLs.

### Required Inputs

* ``--qtlenrich_file1``: QTLEnrich results output file used as input into all types of plotting features (QTLEnrich_{qtl_type}_{trait_name}_{exp_label}_{date}.tsv).

This module requires at least one of the following plotting options:

* ``--box_plot``: Generates box plot.

* ``--violin_plot``: Generates violin plot.

### Optional Inputs

If user is comparing two outputs from **QTLEnrich**:

* ``--qtlenrich_file2``: Second QTLEnrich output file if comparing results from two QTL outputs (e.g., eQTL and sQTL) (QTLEnrich_{qtl_type}_{trait_name}_{exp_label}_{date}.tsv).

* ``--paired_plot``: Indicates user wants to generate a paired box plot or violin plot (e.g., two QTL types).

### Other options:

* --bar_box_violin_option: Options include: 'pi1\_trait', which plots number of trait associations based on an empirical Pi1; 'pi1', which plots the empirical Pi1 rate; 'adjusted\_fold', which plots adjusted fold-enrichment (default: pi1\_trait)

* --label_file1: Label for first QTLEnrich output file (e.g., eQTL). If label is more than one word, e.g. Type 2 Diabetes, please add quotes to variable or string in shell script (default: "QTL type 1"):

	./forest_boxplot_module.R --label_file1 "Type 2 Diabetes"        

* --label_file2: Label for second QTLEnrich output file (e.g., sQTL). For format of input, please see \-\-label\_file1 (default: "QTL type 2")

* --significant_option: Flag specifies generating a plot including only the significant tissues. The significance cutoff is defined in following option:

        * --identify_significant_cutoff: Enrichment p-value cutoff for significant tissues (e.g. Enrichemnt P<0.001). If not indicated, default is Bonferroni cutoff (default: 0.05/number of tissues in QTLEnrich output file).

* --identify_number_significant_tissues: If user wants to plot a susbset of top ranked tissues, ranked based on their adjusted fold-enrichment, e.g. first 10 tissues.

* --remove_subset_tissues: To exclude specific tissues, include text file with list of tissues to exclude. This file is comprised of one column with header 'Tissue (QTL)' and tissue names should match input QTL files. Official GTEx tissue names are listed in: "../data/GTEx\_Tissue\_Names.txt".

        Tissue (QTL)
        Whole_Blood
        Nerve_Tibial
        Brain_Spinal_cord_cervical_c-1

* --gtex_tissue_names: If QTLs are from GTEx, this option uses GTEx tissue names as tissue labels. Tissue names can be found in this file downloaded from the GTEx portal for v8: "../data/gtex_tissue_colors_v8.txt."

* --output_directory: Place plots in directory of choosing (default: Same directory name as that specified for the QTLEnrich results files, e.g., Output\_QTLEnrich\_exp\_label\_DATE).

* ``--exp_label``: Add label to filenames (default: no label).
