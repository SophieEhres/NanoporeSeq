#!/bin/bash
#SBATCH --job-name=restrander
#SBATCH --output=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/restrander/restrander_%j.out
#SBATCH --error=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/restrander/restrander_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --account=def-sauvagm

source ~/.bash_profile

nano_dir="$SCRATCH/nanopore/"
fastq_dir="${nano_dir}/combined_fastq"
stranded_dir="${nano_dir}/restrander"

mkdir -p ${stranded_dir}

config_dir="$HOME/tools/restrander/config"


restrander ${fastq_dir}/${1} \
    ${stranded_dir}/${1%.fastq.gz}_restranded.fastq.gz \
    ${config_dir}/PCB111.json \
        > ${stranded_dir}/${1%.fastq.gz}_stats.json
