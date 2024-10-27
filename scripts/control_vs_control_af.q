#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=af-%a-case2
#SBATCH --output=af-%a-case2.out
#SBATCH --array=16-20
#SBATCH --mem=125G
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP
# make directory for different pvalue thresholds:
     
# generate data for replicates
base_dir='/hpc/group/ochoalab/tt207/platform_bias'
input_value=$((SLURM_ARRAY_TASK_ID))
case=2
mkdir $base_dir/sim$case/rep$input_value/af-test
cd $base_dir/sim$case/rep$input_value/af-test

# calculate AF counts for each subpop/ancestry
file_check="bias_S1.acount"

if [ -f "$file_check" ]; then
# If the file exists, skip the chunk
echo "File $file_check exists. Skipping plink AF calculation"
else
  # If the file does not exist, run the script chunk
  echo "File $file_check does not exist. Running plink AF"
  module load Plink/2.00a2LM
  # # AF for p1
  time plink2 --bfile ../30G_6000n_500000m_controls_bias --keep-fam $base_dir/S1.txt --remove ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p1_S1
  time plink2 --bfile ../30G_6000n_500000m_controls_bias --keep-fam $base_dir/S2.txt --remove ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p1_S2
  time plink2 --bfile ../30G_6000n_500000m_controls_bias --keep-fam $base_dir/S3.txt --remove ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p1_S3
  # AF for p2
  time plink2 --bfile ../30G_6000n_500000m_controls_bias --keep-fam $base_dir/S1.txt --keep ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p2_S1
  time plink2 --bfile ../30G_6000n_500000m_controls_bias --keep-fam $base_dir/S2.txt --keep ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p2_S2
  time plink2 --bfile ../30G_6000n_500000m_controls_bias --keep-fam $base_dir/S3.txt --keep ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p2_S3
  # # 
  # 
  # AF for p1
  # time plink2 --bfile ../30G_6000n_500000m_controls --keep-fam $base_dir/S1.txt --remove ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p1_S1_orig
  # time plink2 --bfile ../30G_6000n_500000m_controls --keep-fam $base_dir/S2.txt --remove ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p1_S2_orig
  # time plink2 --bfile ../30G_6000n_500000m_controls --keep-fam $base_dir/S3.txt --remove ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p1_S3_orig
  # # AF for p2
  # time plink2 --bfile ../30G_6000n_500000m_controls --keep-fam $base_dir/S1.txt --keep ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p2_S1_orig
  # time plink2 --bfile ../30G_6000n_500000m_controls --keep-fam $base_dir/S2.txt --keep ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p2_S2_orig
  # time plink2 --bfile ../30G_6000n_500000m_controls --keep-fam $base_dir/S3.txt --keep ../id_p2.txt --allow-no-sex --nonfounders --freq counts --out p2_S3_orig
  # 
  
  module unload Plink/2.00a2LM

fi

# load R
export LD_LIBRARY_PATH=/opt/apps/rhel8/lapack/lib64:$LD_LIBRARY_PATH
module load R/4.1.1-rhel8
# run AF test
echo "Run AF test to get p-values, OR, and CIs"
time Rscript $base_dir/control_vs_control_af_new.R -r $input_value -s $case -p "1e-02,1e-04,1e-06,1e-08,1e-10,1e-12"
echo "Finished AF test"
module unload R/4.1.1-rhel8
