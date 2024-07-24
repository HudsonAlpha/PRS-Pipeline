# PRS-Pipeline

## Description
This Bash script generates and analyzes three sets of Polygenic Risk Scores (PRS) using the PLINK, PRSice-2, and lassosum software for each individual in the bfiles. It runs the input files through each software to calculate PRS, then analyzes the scores against the phenotypes designated in corresponding fam files. For each fam file, Welch’s t-tests comparing means of PRS for cases and controls are conducted and a corresponding density plot is generated for each PRS program.

## Requirements
- PLINK 1.90
- PRSice-2 (version 2.3.5)
- lassosum (version 0.4.5)

## Installation
Change hard coded paths to PLINK, PRSice-2, lassosum executables to your correct paths. 

GWAS tables must have the following information:
- variant
- allele
- pValue
- beta
- chromosome 
- position

If you are using a different GWAS table, you may need to change the numbers for PLINK’s score function to select the correct columns in this order:
1. Variant id
2. Effect allele
3. Effect size

Additionally, you may need to change the column names following some PRSice-2 flags (snp, a1, pvalue, stat) and lassosum flags (chr, pos, A1, pval, beta) to your corresponding GWAS column names.

## Usage
To run this script, you must provide it with the following four arguments: GWAS table in tsv format, the stem of the PLINK bfiles, a text file including the path to desired fam files separated by a new line, and the covariate file.

```bash
./pd_prs_pipeline.sh <gwas_table>.tsv <bfile_stem> <fam_file_list>.txt <covariate_file>.cov
```

## Arguments
- <gwas_table>.tsv - The GWAS table
- <bfile_stem> - Stem to PLINK executables (.bim & .bed)
- <fam_file_list>.txt - The list of paths to the desired fam files for analysis (the first file will be used for all the PRS calculations)
- <covariate_file>.cov - The covariate file

## Output
- A tsv file for each fam file containing 5 columns: id, PLINK prs, PRSice-2 prs, lassosum prs, and phenotype
- Three density plots (1 for each PRS software) depicting the PRS density by phenotype for each fam file
- PLINK PRS output files (.log, .nopred, .profile)
- PRSice-2 output files (.best, .log, .prsice, .summary, barplot, high-res plot)
- lassosum output files (.lassosum.pipeline.rds, .splitvalidate.rds, .splitvalidate.results.txt, .validate.rds, .validate.results.txt)

## References
