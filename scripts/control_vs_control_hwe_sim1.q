#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=hwe-%a-case1
#SBATCH --output=hwe-%a-case1.out
#SBATCH --array=1
#SBATCH --mem=125G
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP
# make directory for different pvalue thresholds:
     
# generate data for replicates
base_dir='/hpc/group/ochoalab/tt207/platform_bias'
input_value=$((SLURM_ARRAY_TASK_ID))
case=1
PVAL=1e-10

mkdir $base_dir/sim$case/rep$input_value/hwe-test
cd $base_dir/sim$case/rep$input_value/hwe-test
mkdir $base_dir/sim$case/rep$input_value/hwe-test/$PVAL
cd $PVAL

# # get hardy pvals
module load Plink/2.00a2LM
echo 'calculate hwe pvalues'
plink2 --bfile ../../1G_6000n_500000m_controls --hardy --nonfounders --out bias_hwe
module unload Plink/2.00a2LM
echo 'classify sig snps into remove and flip'
export LD_LIBRARY_PATH=/opt/apps/rhel8/lapack/lib64:$LD_LIBRARY_PATH
module load R/4.1.1-rhel8
time Rscript $base_dir/hardy_step1.R -r $input_value -s $case -p "1e-10"
#
module load Plink/1.90
echo "flip selected SNPs"
time plink --bfile ../../1G_6000n_500000m_controls \
    --flip hwe_flip_fwd.txt \
    --flip-subset ../../id_p2.txt \
    --extract hwe_flip_fwd.txt \
    --make-bed --out hwe_flip_reverse
module unload Plink/1.90
echo "calculate hwe to flips in reverse"
module load Plink/2.00a2LM
plink2 --bfile hwe_flip_reverse --hardy --nonfounders --out hwe_flip_reverse
module unload Plink/2.00a2LM
#
echo "compare fwd and rev flip p-values"
time Rscript $base_dir/hardy_step2.R -r $input_value -s $case -p "1e-10"

echo "perform finalized flip and remove"
module load Plink/1.90
time plink --bfile ../../1G_6000n_500000m_controls --flip hwe_flip_final.txt --flip-subset ../../id_p2.txt \
    --exclude hwe_remove_final.txt \
    --make-bed --out hwe_clean
module unload Plink/1.90

echo "calculate AF of new files"
# # # calculate AF - clean
module load Plink/2.00a2LM
time plink2 --bfile hwe_clean --keep ../../id_p1.txt --keep-fam $base_dir/S1.txt --freq counts --nonfounders --out p1_S1
time plink2 --bfile hwe_clean --keep ../../id_p1.txt --keep-fam $base_dir/S2.txt --freq counts --nonfounders --out p1_S2
time plink2 --bfile hwe_clean --keep ../../id_p1.txt --keep-fam $base_dir/S3.txt --freq counts --nonfounders --out p1_S3
# AF for tgp
time plink2 --bfile hwe_clean --keep ../../id_p2.txt --keep-fam $base_dir/S1.txt --freq counts --nonfounders  --out p2_S1
time plink2 --bfile hwe_clean --keep ../../id_p2.txt --keep-fam $base_dir/S2.txt --freq counts --nonfounders  --out p2_S2
time plink2 --bfile hwe_clean --keep ../../id_p2.txt --keep-fam $base_dir/S3.txt --freq counts --nonfounders  --out p2_S3
module unload Plink/2.00a2LM



