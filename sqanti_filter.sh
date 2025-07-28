#!/bin/bash
#SBATCH --job-name=sqantiFilter
#SBATCH --output=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/sqantiFilter/sqantiFilter_%j.out
#SBATCH --error=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/sqantiFilter/sqantiFilter_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --account=def-sauvagm

module load python/3.11.5 gcc/12.3 minimap2/2.28 scipy-stack/2025a samtools/1.20 star/2.7.11b r/4.4.0

sqanti_dir="/home/p1044860/tools/SQANTI3"


source ~/.bash_profile
source "${sqanti_dir}/sqanti3_env/bin/activate"

#Set up directories
nano_dir="${SCRATCH}/nanopore"
anno_sqanti_dir="${nano_dir}/sqantiQC/minimap2_restrander"
genome_dir="/home/p1044860/projects/def-sauvagm/ref_genomes/human/CHM13"
out_dir="${nano_dir}/sqantiFilter/minimap2_restrander"
sam_dir="${nano_dir}/minimap2_restrander"

mkdir -p ${out_dir}

#Set up target file

sample_dir="${1}"
name="${1%.sam}"
classification_file="${sample_dir}_classification.txt"
gtf_file="${sample_dir}_corrected.gtf"
sam_file=$(ls ${sam_dir} | grep -e ".sam" | grep -e "${name}" | grep -v ".bam")

## run sqanti
python ${sqanti_dir}/sqanti3_filter.py \
    rules \
    --sqanti_class ${anno_sqanti_dir}/${sample_dir}/${classification_file} \
    --filter_sam ${sam_dir}/${sam_file} \
    --filter_gtf ${anno_sqanti_dir}/${sample_dir}/${gtf_file} \
    -c 4 \
    -o ${name} \
    -d ${out_dir}/${name}
