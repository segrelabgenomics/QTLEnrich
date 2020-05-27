# README of Examples folder

To run **QTLEnrich**, one needs to input a GWAS summary statistics results file and a set of QTLs, both of which *must* be in the same genome build, e.g. GRCh38.  We provide an example run of **QTLEnrich** using a sample Height GWAS in genome build GRCh37 and whole blood eQTLs and sQTLs from GTEx release v8.

# Steps to obtain sample GWAS and GTEx eQTL and sQTL datasets and liftOver GWAS to build GRCh38.

All example input files should be saved in the current **examples/** folder.

## 1.  Download sample GWAS

We downloaded the Height GWAS computed by the Neale lab at MGH and Broad Institute, using the UK Biobank data, as described here: http://www.nealelab.is/uk-biobank.
The code used to generate this height GWAS is publicly available: https://github.com/Nealelab/UK_Biobank_GWAS.

```
wget https://www.dropbox.com/s/ou12jm89v74k55e/50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz?dl=0 -O 50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz
zcat 50_irnt.gwas.imputed_v3.both_sexes.tsv.bgz > 50_irnt.gwas.imputed_v3.both_sexes.tsv 
```

The GWAS file has to have at least the chr number and position of each variant and the association p-value. For more information, see README within *src* directory.  

## 2. liftOver GWAS from genome build 37 to build 38 and process columns and headers to conform with the **QTLEnrich** input GWAS file requirements (see **src/README.md**)

We provide a wrapper script to liftover the sample UK Biobank Height GWAS from genome build GRCh37 to build GRCh38. This script was written in Python 3.

```
./run_liftOver_gwas.sh
```

This script utilizes the liftOver tool and chain file from UCSC that can be downloaded from:

liftOver tool: [http://hgdownload.soe.ucsc.edu/admin/exe/](http://hgdownload.soe.ucsc.edu/admin/exe/)

liftOver chain file from genome build GRCh37 to GRCh38, hg19ToHg38.over.chain.gz: [http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/)

**Note:** This liftOver wrapper scrtipt was writted for the UK Biobank GWAS output files from the Neale lab, but can be adapted for other GWAS formats. 

## 3. Download QTL files

GTEx v8 eQTLs and sQTLs can be downloaded from the [GTEx portal](https://www.gtexportal.org/home/datasets) or from [Github](https://github.com/broadinstitute/gtex-v8):

### GTEx v8 eQTLs
```
wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL.tar && tar xf GTEx_Analysis_v8_eQTL.tar && rm GTEx_Analysis_v8_eQTL.tar
```

### GTEx v8 sQTLs
```
wget https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_sQTL.tar && tar xf GTEx_Analysis_v8_sQTL.tar && rm GTEx_Analysis_v8_sQTL.tar
```

## Sample eQTL and sQTL files used here:

* Whole_Blood.v8.egenes.txt.gz - best eQTL per gene in whole blood from GTEx release v8
* Whole_Blood.v8.sgenes.txt.gz - best sQTL per gene in whole blood from GTEx release v8

## 4. Run **QLTEnrich** using the command in *src* folder:
```
./sample_qtlenrich_command.sh
```

A detailed description of the inputs and outputs of **QTLEnrich** is given in **src/README.md**

## 5. Output files from sample run that we performed can be found in subdirectory: 

Output_QTLEnrich_example_Mar02_2020/

