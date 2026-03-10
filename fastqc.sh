#!/bin/bash

##SBATCH --account=def-sauvagm
##SBATCH --mem-per-cpu=20G
##SBATCH -n 10
##SBATCH --time=20:00

fastqdir="/home/p1044860/projects/def-sauvagm/data/nanopore_seq/2024_04_11_MCF10a_TGFb_timecourse/fastq"
outdir="/home/p1044860/projects/def-sauvagm/data/nanopore_seq/2024_04_11_MCF10a_TGFb_timecourse/fastqc"

file=$1 ## Get from bash input

module load fastqc



fastqc --outdir ${outdir} \
	--nano \
	-t 10 \
	--memory 1000 \
	${fastqdir}/${file} 
