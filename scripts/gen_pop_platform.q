#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=genpop_case2
#SBATCH --output=genpop_case2_%a.out
#SBATCH --array=16-20
#SBATCH --mem=120G
#SBATCH --ntasks-per-node=60
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# simulation scenario
case='2'
base_dir='/hpc/group/ochoalab/tt207/platform_bias'
# generate data for 10 replicates
input_value=$((SLURM_ARRAY_TASK_ID))
filename='1G_6000n_500000m_controls'

export LD_LIBRARY_PATH=/opt/apps/rhel8/lapack/lib64:$LD_LIBRARY_PATH
module load R/4.1.1-rhel8
# echo "create plink files for 2 platforms"
# # simulate data from 6000 individuals, then later divide into two platforms
# 
# ##### SIM1 #####
# time Rscript gen_pop_sim1.R -n $input_value -s 1
# echo "finishing running gen_pop.R"
# 
# time Rscript bias_geno.R -n $input_value -s $case -f 1G_6000n_500000m_controls
# echo "finishing running bias_geno.R"

# create link file names for bim and fam file
# rm ./sim$case/rep$input_value/1G_6000n_500000m_controls.bed
# mv ./sim$case/rep$input_value/1G_6000n_500000m_controls_bias.bed ./sim$case/rep$input_value/1G_6000n_500000m_controls.bed
################

##### SIM2 ######
#split sim1 data by subpop for sim2
mkdir sim1/rep$input_value/subpop
module load Plink/1.90
time plink --keep-allele-order --bfile ./sim1/rep$input_value/${filename} --keep-fam S1.txt --make-bed --out ./sim1/rep$input_value/subpop/1G_2000n_500000m_S1
time plink --keep-allele-order --bfile ./sim1/rep$input_value/${filename} --keep-fam S2.txt --make-bed --out ./sim1/rep$input_value/subpop/1G_2000n_500000m_S2
time plink --keep-allele-order --bfile ./sim1/rep$input_value/${filename} --keep-fam S3.txt --make-bed --out ./sim1/rep$input_value/subpop/1G_2000n_500000m_S3
time Rscript gen_pop_sim2_step1.R -n $input_value -s 2
echo "finished step1"
time Rscript gen_pop_sim2_step2.R -n $input_value -s 2
echo "finished step2"
module load Plink/1.90
plink --bfile ./sim2/rep$input_value/subpop/30G_2000n_500000m_S1 --bmerge ./sim2/rep$input_value/subpop/30G_2000n_500000m_S2 --make-bed --out ./sim2/rep$input_value/merged_temp
plink --bfile ./sim2/rep$input_value/merged_temp --bmerge ./sim2/rep$input_value/subpop/30G_2000n_500000m_S3 --make-bed --out ./sim2/rep$input_value/30G_6000n_500000m_controls
module unload Plink/1.90
rm ./sim2/rep$input_value/merged_temp*
time Rscript gen_pop_sim2_step3.R -n $input_value -s 2
echo "finished step3"

time Rscript bias_geno.R -n $input_value -s 2 -f 30G_6000n_500000m_controls
echo "finishing running bias_geno.R"
rm -r -d ./sim2/rep$input_value/subpop
mkdir ./sim2/rep$input_value/saige
mkdir ./sim2/rep$input_value/af-test
# rm ./sim$case/rep$input_value/30G_6000n_500000m_controls.*
# for file in ./sim$case/rep$input_value/30G_6000n_500000m_controls_bias*; do
#   new_file=$(echo "$file" | sed 's/_bias//')
#   mv "$file" "$new_file"
# done
# 
echo "remove subpop directory and create saige and af-test dir"
echo "remove original file and rename bias file"