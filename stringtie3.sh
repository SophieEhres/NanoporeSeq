#!/bin/bash
#SBATCH --job-name=stringtie
#SBATCH --output=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/stringtie/stringtie_%j.out
#SBATCH --error=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/stringtie/stringtie_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8
#SBATCH --account=def-sauvagm

module load samtools

# Set up directories and files
nano_dir="${SCRATCH}/nanopore"
bam_dir="${nano_dir}/minimap2_restrander"
SR_bam_dir="${SCRATCH}/MCF10a_illumina/bam"
out_dir="${nano_dir}/stringtie_3/minimap2_restrander"

file="${1}"
LR_bam="${bam_dir}/${file}"


## get values to look for SR file
time=$(echo ${file} | cut -d "_" -f1)
condition=$(echo ${file} | cut -d "_" -f2)
replicate=$(echo ${file} | cut -d "_" -f3 | cut -d "." -f1)


case "$time" in ### get corresponding time value
    0h) SR_time="D0" ;;
    12h) SR_time="12hrs" ;;
    24h) SR_time="D1" ;;
    48h) SR_time="D2" ;;
    96h) SR_time="D4" ;;
esac

case "$condition" in ### get corresponding condition value
    Non-Treated) SR_condition="NT" ;;
    TGFB1) SR_condition="TGF" ;;
esac

case "$replicate" in ### get corresponding replicate value
    rep1) SR_rep="1" ;;
    rep2) SR_rep="2" ;;
    rep3) SR_rep="3" ;;
esac

## get SR file name
if [[ $time == "0h" ]]; then
    SR_file_name="${SR_time}_${SR_rep}"
else
    SR_file_name="${SR_time}_${SR_condition}_${SR_rep}"
fi

# find SR file
SR_file=$(ls ${SR_bam_dir} | grep -e ${SR_file_name} | grep -v "bai" | grep -v "tmp" | head -n 1)
SR_bam="${SR_bam_dir}/${SR_file}"

echo -e "${file} ......... ${SR_file}"

if [[ "$file" == *_sorted.bam ]]; then
    exit
else
    samtools sort -@ 8 -o ${LR_bam%.bam}_sorted.bam ${LR_bam}
    LR_bam="${LR_bam%.bam}_sorted.bam"
fi

### files needed for stringtie

ref_gtf="/home/p1044860/projects/def-sauvagm/ref_genomes/human/CHM13/Homo_sapiens-GCA_009914755.4-2022_07-genes_chr.gtf"

echo "Using SR BAM: $SR_bam"
echo "Using LR BAM: $LR_bam"
echo "Saving output to: ${out_dir}/${file%.bam}.gtf"


stringtie ${SR_bam} ${LR_bam} \
    --mix \
    -p 8 \
    -G ${ref_gtf} \
    -l MCF10a_mix_minimapOnly \
    -o ${out_dir}/${file%.bam}.gtf 
