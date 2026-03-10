#!/bin/bash
#SBATCH --account=def-sauvagm
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --output=/home/p1044860/projects/def-sauvagm/data/nanopore_seq/2024_04_11_MCF10a_TGFb_timecourse/scripts/bwtool_out/bwtool_%j.out
#SBATCH --error=/home/p1044860/projects/def-sauvagm/data/nanopore_seq/2024_04_11_MCF10a_TGFb_timecourse/scripts/bwtool_out/bwtool_%j.err

module load kentutils bedops


gff="/home/p1044860/projects/def-sauvagm/data/nanopore_seq/2024_04_11_MCF10a_TGFb_timecourse/bwtool/SE_MCF10a_Nanopore_Dorado_Restrander_Minimap2_Stringtie3_Sqanti_merge_hg38_liftover.gtf"
bed="${gff%.gtf}.bed"
format_bed="${bed%.bed}_formatted.bed"
bw="/home/p1044860/projects/def-sauvagm/ref_genomes/human/hg38/gerp_conservation_scores.homo_sapiens.GRCh38.bw"

bwdir="/home/p1044860/projects/def-sauvagm/data/nanopore_seq/2024_04_11_MCF10a_TGFb_timecourse/bwtool"

mkdir -p ${bwdir}

gff2bed < ${gff} > ${bed}
cat ${bed} | awk '{gsub(/^chr/, "", $1); if($5 == ".") $5 = 0; print $1,$2,$3,$4,$5,$6}' | tr ' ' '\t' > ${format_bed}
bwtool summary -header ${format_bed} ${bw} ${bwdir}/SE_MCF10a_Nanopore_Dorado_Restrander_Minimap2_Stringtie3_Sqanti_merge_hg38_liftover_gerpscore.txt
