#!/bin/bash -l
#SBATCH --job-name=T2D_PRS
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --cpus-per-task=1

### Example batch script for PRSice. All paths need updating. All file names should
### be checked.

module add r


    /path_to_PRSice_software/PRSice/PRSice_linux \
    --a1 A1 \
    --a2 A2 \
    --bar-levels 5e-08,1e-05,0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5,1 \
    --base /path_to_disc_t2d_sumstats.gz \
    --binary-target T \
    --clump-kb 250kb \
    --clump-p 1.000000 \
    --clump-r2 0.100000 \
    --cov /path_to_covariates_to_adjust_for/applicationnumber_covariates.txt \
    --cov-col @PC[1-6],batch,assessment_centre \
    --cov-factor batch,assessment_centre \
    --extract /path_to_postQC_snplist/UKB.postQC.EUR.Autosome.MAF_01.GENO_02.HWE_1e-8.snplist \
    --fastscore  \
    --keep /analysis_dir/data/EUR_postQC.fam \
    --num-auto 22 \
    --out /analysis_dir/PRSice_out \
    --pheno /analysis_dir/data/EURt2d_pheno_file_for_PRSice.txt \
    --pvalue P \
    --seed 527288700 \
    --snp SNP \
    --stat BETA \
    --target /path_to_your_ukbiobank_genotyped_data/ukbprj_binary_pre_qc,/path_to_your_ukbiobank_genotyped_data/ukbprj.fam \
    --thread 1
