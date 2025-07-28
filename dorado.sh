#!/bin/bash
#SBATCH --job-name=dorado
#SBATCH --output=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/dorado_%j.out
#SBATCH --error=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/dorado_%j.err
#SBATCH --time=3-00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=1
#SBATCH --account=def-sauvagm

module load StdEnv/2023 dorado/0.8.3 

nano_dir="/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse"
pod5_root="${nano_dir}/pod5"
sample_sheet="${nano_dir}/scripts/sample_sheet_MSauvageau_SQuaggin_exp.csv"
out_dir="$SCRATCH/nanopore/dorado"

# Directory name passed from command line
dir_name="$1"
pod5_sample_dir="${pod5_root}/${dir_name}"

mkdir -p "${out_dir}/${dir_name}"

dorado basecaller \
    /home/p1044860/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
    "${pod5_sample_dir}" \
    --kit-name SQK-NBD114-24 \
    --no-trim \
    --sample-sheet "${sample_sheet}" \
    --verbose \
    --estimate-poly-a \
    --min-qscore 7 \
    --emit-sam \
    > "${out_dir}/${dir_name}.bam"

dorado summary "${out_dir}/${dir_name}.bam" > "${out_dir}/summary_${dir_name}.tsv"

dorado demux \
    "${out_dir}/${dir_name}.bam" \
    --kit-name SQK-NBD114-24 \
    --sample-sheet "${sample_sheet}" \
    -t 16 \
    -v \
    --emit-fastq \
    --no-trim \
    --emit-summary \
    -o "${out_dir}"

