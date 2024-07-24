#!/bin/bash

#SBATCH -c 8
#SBATCH --mem=64G

# Message to print if # of inputs do not match required
if [ $# -ne 4 ];
then
    echo "Usage: <gwas_sumstats> <bfile> <famfile_list> <covariatefile>"
    exit 1
fi

# Load necessary modules
module load cluster/plink/1.90
source /cluster/home/dkim/y/etc/profile.d/micromamba.sh
micromamba activate /cluster/home/dkim/y/envs/park_prs

# Assign variables
gwas_ss=$1
bfile=$2
fam_list=$3
cov=$4

fam=$(head -n 1 $fam_list)

prsice_script="/cluster/home/dkim/PRSice_2.3.5/PRSice.R"
lassosum_script="/cluster/home/dkim/y/envs/park_prs/lib/R/library/lassosum/lassosum"

working_dir=$(pwd)

output_basename=$(basename "$bfile")

# Generate PLINK prs
plink --bfile $bfile \
    --fam $fam \
    --score $gwas_ss 1 2 4 header no-mean-imputation \
    --out ${output_basename}_plink \
    --memory 64000

# Generate PRSice-2 prs
Rscript $prsice_script \
    --prsice /cluster/home/dkim/PRSice_2.3.5/PRSice_linux \
    --base $gwas_ss \
    --target $bfile,$fam \
    --snp variant \
    --a1 allele \
    --pvalue pValue \
    --stat beta \
    --binary-target T \
    --beta \
    --thread 8 \
    --cov $cov \
    --out ${output_basename}_prsice

# Generate lassosum prs
echo -e "FID\tIID\tPhenotype" > lasso_pheno.file
cat $fam | awk '{ print $1 "\t" $2 "\t" $6 }' >> lasso_pheno.file

$lassosum_script \
    --data $gwas_ss \
    --chr chromosome \
    --pos position \
    --A1 allele \
    --pval pValue \
    --beta beta \
    --n 482730 \
    --test.bfile $bfile \
    --LDblocks EUR.hg38 \
    --pheno lasso_pheno.file \
    --covar $cov \
    --nthreads 8 \
    --out ${output_basename}_lassosum

# remove pheno file
rm lasso_pheno.file

# Run analysis and produce plots based on results
Rscript /cluster/home/dkim/prs_analysis_all_fam.r $fam_list ${output_basename}_plink.profile ${output_basename}_prsice.best ${output_basename}_lassosum.splitvalidate.results.txt $working_dir
