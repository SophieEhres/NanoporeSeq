#!/bin/bash
#SBATCH --job-name=dorado
#SBATCH --output=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/minimap2/minimap2_%j.out
#SBATCH --error=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/minimap2/minimap2_%j.err
#SBATCH --time=3:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=16
#SBATCH --account=def-sauvagm

nano_dir="$SCRATCH/nanopore"
fastq_dir="${nano_dir}/restrander"
bam_dir="${nano_dir}/minimap2_restrander"
genome_dir="/home/p1044860/projects/def-sauvagm/ref_genomes/human/CHM13/fasta"
genome="${genome_dir}/chm13v2.0.fa"
file=${1}

module load minimap2

mkdir -p ${bam_dir}

minimap2 -ax splice \
    -t 16 \
    -o ${bam_dir}/${file%.fastq.gz}.sam \
    ${genome} \
    ${fastq_dir}/${file}
    




