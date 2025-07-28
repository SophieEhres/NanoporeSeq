#!/bin/bash
#SBATCH --job-name=pychopper_reads
#SBATCH --output=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/pychopper/pychopper_%j.out
#SBATCH --error=/home/p1044860/projects/def-sauvagm/data/nanopore-seq/2024_SE_MCF10a_EMT_Timecourse/scripts/pychopper/pychopper_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=12
#SBATCH --account=def-sauvagm


source /home/p1044860/tools/pychopper/pychopper_env/bin/activate
module load python/3.11.5 parasail hmmer

nano_dir="$SCRATCH/nanopore"
fastq_dir="${nano_dir}/combined_fastq"
pychopper_dir="${nano_dir}/pychopper"

mkdir -p ${pychopper_dir}

name=${1%.fastq.gz}

pychopper -k LSK114 \
    -r ${pychopper_dir}/${name}_report_pHMM.pdf \
    -u ${pychopper_dir}/${name}_unclassified_pHMM.fq \
    -w ${pychopper_dir}/${name}_rescued_pHMM.fq \
    -t 12 \
    ${fastq_dir}/${1} \
    ${pychopper_dir}/${name}_full_length_pHMM.fq

cat ${pychopper_dir}/${name}_full_length_pHMM.fq ${pychopper_dir}/${name}_rescued_pHMM.fq > ${pychopper_dir}/${name}_concat_pHMM.fq



pychopper -k LSK114 \
    -r ${pychopper_dir}/${name}_report_edlib.pdf \
    -u ${pychopper_dir}/${name}_unclassified_edlib.fq \
    -w ${pychopper_dir}/${name}_rescued_edlib.fq \
    -t 12 \
    ${fastq_dir}/${1} \
    ${pychopper_dir}/${name}_full_length_edlib.fq

cat ${pychopper_dir}/${name}_full_length_edlib.fq ${pychopper_dir}/${name}_rescued_edlib.fq > ${pychopper_dir}/${name}_concat_edlib.fq
