#!/bin/bash
#SBATCH --job-name=stringtie_merge
#SBATCH --output=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/stringtie_merge/stringtie_merge_%j.out
#SBATCH --error=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/stringtie_merge/stringtie_merge_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --account=def-sauvagm

module load samtools

# set up directories
nano_dir="${SCRATCH}/nanopore"
rescue_dir="${nano_dir}/sqantiRescue"
genome_dir="/home/p1044860/projects/def-sauvagm/ref_genomes/human/CHM13"

out_dir="${nano_dir}/stringtie_merge_from_rescue"

mkdir -p ${out_dir}

#set up required files

ref_gtf="${genome_dir}/Homo_sapiens-GCA_009914755.4-2022_07-genes_chr_filterFrames.gtf" ## Use the one with transcripts were removed if they had same start/stop codon [no frame error gffread]

gtf_files=$(ls ${rescue_dir}/minimap2_restrander/*/*.gtf | grep -e "reorder" | tr '\n' ' ') ### reordered gtf to correct negative strand start/end

echo -e "gtf files \n ${gtf_files}"

stringtie --merge \
    -G ${ref_gtf} \
    -o ${out_dir}/SE_MCF10a_Nanopore_Dorado_Restrander_Minimap2_Stringtie3_Sqanti_merge.gtf \
    -l MCF10a \
    ${gtf_files}
