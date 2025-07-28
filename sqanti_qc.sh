#!/bin/bash
#SBATCH --job-name=sqantiQC
#SBATCH --output=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/sqantiQC/sqantiQC_%j.out
#SBATCH --error=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/sqantiQC/sqantiQC_%j.err
#SBATCH --time=00:30:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --account=def-sauvagm

module load python/3.11.5 gcc/12.3 minimap2/2.28 scipy-stack/2025a samtools/1.20 star/2.7.11b r/4.4.0

sqanti_dir="/home/p1044860/tools/SQANTI3"


source ~/.bash_profile
source "${sqanti_dir}/sqanti3_env/bin/activate"

#Set up directories
nano_dir="${SCRATCH}/nanopore"
anno_dir="${nano_dir}/stringtie_3/minimap2_restrander_c9"
genome_dir="/home/p1044860/projects/def-sauvagm/ref_genomes/human/CHM13"
out_dir="${nano_dir}/sqantiQC/minimap2_restrander_c9"

mkdir -p ${out_dir}

#Set up required files
ref_gtf="${genome_dir}/Homo_sapiens-GCA_009914755.4-2022_07-genes_chr_filterFrames.gtf" ## Use the one with transcripts were removed if they had same start/stop codon [no frame error gffread]
ref_genome="${genome_dir}/fasta/chm13v2.0.fa"
polyA="${sqanti_dir}/data/polyA_motifs/mouse_and_human.polyA_motif.txt"

#Set up target file

file="${1}"
name="${file%.gtf}"

## run sqanti
python ${sqanti_dir}/sqanti3_qc.py \
    --isoforms "${anno_dir}"/"${file}" \
    --refGTF "${ref_gtf}" \
    --refFasta "${ref_genome}" \
    --polyA_motif_list ${polyA} \
    -t 4 \
    -o "${name}" \
    -d "${out_dir}"/${name} \
    --report both 

