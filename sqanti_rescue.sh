#!/bin/bash
#SBATCH --job-name=sqantiRescue
#SBATCH --output=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/sqantiRescue/sqantiRescue_%j.out
#SBATCH --error=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/sqantiRescue/sqantiRescue_%j.err
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
qc_dir="${nano_dir}/sqantiQC/minimap2_restrander"
filter_dir="${nano_dir}/sqantiFilter/minimap2_restrander"
genome_dir="/home/p1044860/projects/def-sauvagm/ref_genomes/human/CHM13"
out_dir="${nano_dir}/sqantiRescue/minimap2_restrander"

mkdir -p ${out_dir}

name=${1}

#Set up required files
ref_gtf="${genome_dir}/Homo_sapiens-GCA_009914755.4-2022_07-genes_chr_filterFrames.gtf" ## Use the one with transcripts were removed if they had same start/stop codon [no frame error gffread]
ref_genome="${genome_dir}/fasta/chm13v2.0.fa"

## run sqanti
python ${sqanti_dir}/sqanti3_rescue.py \
    rules \
    --filter_class ${filter_dir}/${name}/*_classification.txt \
    --refGTF ${ref_gtf} \
    --refFasta ${ref_genome} \
    --rescue_gtf ${filter_dir}/${name}/*.filtered.gtf \
    -o ${name} \
    -d ${out_dir}/${name} \
    -c 4

