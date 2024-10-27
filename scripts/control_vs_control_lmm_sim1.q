#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=saige_1e-02-sim1
#SBATCH --output=rep%a_1e-02_sim1.out
#SBATCH --array=1-10
#SBATCH --mem=125G
#SBATCH --ntasks-per-node=24
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP
# make directory for different pvalue thresholds:

# generate data for replicates
base_dir='/hpc/group/ochoalab/tt207/platform_bias'
input_value=$((SLURM_ARRAY_TASK_ID)) # rep
PVAL=1e-02
case=1
cd $base_dir/sim$case/rep$input_value/saige

module load Plink/1.90

# create pval directory if it doesnt exist
if [ -d "${PVAL}" ]; then
  # if the directory exists, delete
  rm -r "${PVAL}"
fi
#
mkdir ${PVAL}
cd ${PVAL}
# 
# # load R
export LD_LIBRARY_PATH=/opt/apps/rhel8/lapack/lib64:$LD_LIBRARY_PATH
module load R/4.1.1-rhel8
module load Plink/1.90
# 
# # check if iteration 0 saige output exists, if not, run saige for iteration 0 
ITER=0
# # 
file_path="$base_dir/sim$case/rep$input_value/saige/platform_0_output.txt"
if [ -f "$file_path" ]; then
  echo "saige output for iteration 0 already exists. Skipping initial saige execution"
else
  echo "saige running for iteration 0"
  time Rscript $base_dir/saige_step1.R -i $ITER -p $base_dir/sim$case/rep$input_value/saige/${PVAL} -s $case
  time Rscript $base_dir/saige_step2.R -i $ITER -p $base_dir/sim$case/rep$input_value/saige/${PVAL} -s $case
fi

# # # # # # # create summary file 'remove' for each iteration
summary_file="iteration_summary_2.txt"
echo -e "ITER\tremove" > "$summary_file"



while :; do
    suffix=$ITER
    next_suffix=$((ITER + 1))
    time Rscript $base_dir/phase1_remove.R -i $ITER -p $PVAL -r $input_value -s $case

    # run plink to remove/flip snps
    # in PVAL directory

    if [ "$ITER" -eq 0 ]; then

        # Count rows in flip and remove files
        remove_count=$(($(wc -l < platform_snps_phase1.txt) - 1))
        echo -e "$ITER\t$remove_count" >> "$summary_file"
        if [ "$remove_count" -le 0 ]; then
        echo "remove_count is zero. Stopping the loop."
        echo -e "$ITER\t$remove_count" >> "$summary_file"
        break
        fi

        time plink --bfile ../../1G_6000n_500000m_controls \
          --exclude platform_snps_phase1.txt \
          --keep-allele-order --make-bed --out 1G_6000n_500000m_controls_$next_suffix

    else

        remove_count=$(($(wc -l < platform_snps_phase1_$suffix.txt) - 1))

        # Ensure flip_count and remove_count are not negative
        remove_count=$((remove_count < 0 ? 0 : remove_count))
        # Append to the summary file
        echo -e "$ITER\t$remove_count" >> "$summary_file"

        if [ "$remove_count" -le 0 ]; then
        echo "remove_count is less than or equal to zero. Stopping the loop."
        break
        fi

        time plink --bfile 1G_6000n_500000m_controls_$suffix \
          --exclude platform_snps_phase1_$suffix.txt \
          --keep-allele-order --make-bed --out 1G_6000n_500000m_controls_$next_suffix

        fi

        # saige
        time Rscript $base_dir/saige_step1.R -i $next_suffix -p $base_dir/sim$case/rep$input_value/saige/${PVAL} -s $case
        time Rscript $base_dir/saige_step2.R -i $next_suffix -p $base_dir/sim$case/rep$input_value/saige/${PVAL} -s $case

        # Increment ITER for the next iteration
        ITER=$next_suffix
done

# # # 
# SECOND PHASE
echo "start phase 2"
time Rscript $base_dir/phase2_flip.R -p $PVAL -r $input_value -s $case
# 
# # # flip and remove snps from original base plink file
# 
time plink --bfile ../../1G_6000n_500000m_controls \
  --flip platform_phase2_flip.txt \
  --flip-subset $base_dir/p2-control-id.txt \
  --exclude platform_phase2_remove.txt \
  --keep-allele-order --make-bed --out 1G_6000n_500000m_controls_phase2
# 
# # saige
time Rscript $base_dir/saige_step1.R -i phase2 -p $base_dir/sim$case/rep$input_value/saige/${PVAL} -s $case
time Rscript $base_dir/saige_step2.R -i phase2 -p $base_dir/sim$case/rep$input_value/saige/${PVAL} -s $case
# 
# # identify sig snps from phase 2 output
time Rscript $base_dir/phase2_remove.R -i phase2 -p $PVAL -r $input_value -s $case
remove_count=$(($(wc -l < platform_sig_snps_remove_phase3_iter0.txt) - 1))
remove_count=$((remove_count < 0 ? 0 : remove_count))
echo -e "phase2\t$remove_count" >> "$summary_file"
if [ "$remove_count" -le 0 ]; then
echo "remove_count is less than or equal to zero. Stopping the rest of the script."
exit 0
fi

echo "start phase 3" 

ITER=0
while :; do
    suffix=$ITER
    prev_suffix=$((ITER - 1))
    next_suffix=$((ITER + 1))
#
    # run plink iterations to remove
    # in PVAL directory
    if [ "$ITER" -eq 0 ]; then
        echo 'iter 0'
        time plink --bfile 1G_6000n_500000m_controls_phase2 \
        --exclude platform_sig_snps_remove_phase3_iter0.txt \
        --keep-allele-order --make-bed --out 1G_6000n_500000m_controls_phase3_0
    else
        echo 'iter $ITER'
        time plink --bfile 1G_6000n_500000m_controls_phase3_$prev_suffix \
        --exclude platform_sig_snps_remove_phase3_iter$suffix.txt \
        --keep-allele-order --make-bed --out 1G_6000n_500000m_controls_phase3_$suffix
    fi
      
        # Run SAIGE for the next iteration
        time Rscript $base_dir/saige_step1.R -i phase3_$suffix -p $base_dir/sim$case/rep$input_value/saige/${PVAL} -s $case
        time Rscript $base_dir/saige_step2.R -i phase3_$suffix -p $base_dir/sim$case/rep$input_value/saige/${PVAL} -s $case

        # Identify significant SNPs
        time Rscript $base_dir/phase2_remove.R -i $suffix -p $PVAL -r $input_value -s $case

        remove_count=$(($(wc -l < platform_sig_snps_remove_phase3_iter${next_suffix}.txt) - 1))
        remove_count=$((remove_count < 0 ? 0 : remove_count))
        echo "Number of lines (excluding header) in platform_sig_snps_remove_phase3_iter${next_suffix}.txt: $remove_count"
        echo -e "$next_suffix\t$remove_count" >> $summary_file

        if [ "$remove_count" -le 0 ]; then
          echo "remove_count is $remove_count. Stopping the loop."
          break
        fi
        # Increment ITER for the next iteration
        ITER=$next_suffix
done

# remove itermediate files
#rm 1G_6000n_500000m_controls_*

module unload R/4.1.1-rhel8

