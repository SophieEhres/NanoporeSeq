#!/bin/bash
#SBATCH --account=def-sauvagm
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
#SBATCH --time=6:00:00
#SBATCH --output=/home/p1044860/links/projects/def-sauvagm/shared/raw_sequencing_data/nanopore-seq/2024_04_11_MCF10a_TGFb_timecourse/scripts/isoquant_out_gabParam/isoquant_%j.out
#SBATCH --error=/home/p1044860/links/projects/def-sauvagm/shared/raw_sequencing_data/nanopore-seq/2024_04_11_MCF10a_TGFb_timecourse/scripts/isoquant_out_gabParam/isoquant_%j.err

echo "Started at: $(date)"
echo "Node: $(hostname)"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
echo "nproc says: $(nproc)"
echo "TMPDIR is: $SLURM_TMPDIR"

# Load required modules
module load StdEnv/2023 samtools/1.22.1 scipy-stack/2025a minimap2/2.28

# Activate environment
isoquant_dir="/home/p1044860/tools/IsoQuant_5.10/IsoQuant"
source "${isoquant_dir}/isoquant_env_5_10/bin/activate"

# Directories
reference_dir="/home/p1044860/links/projects/def-sauvagm/shared/ref_genomes/human/v47"
nano_dir="/home/p1044860/links/projects/def-sauvagm/shared/raw_sequencing_data/nanopore-seq/2024_04_11_MCF10a_TGFb_timecourse"
bam_dir="${nano_dir}/minimap2_restrander_hg38"
final_out_dir="${SCRATCH}/isoquant_dorado_restranded_minimap2_gab_param"

mkdir -p ${final_out_dir}

file=${1}

# ------------------------------
# Copy input files to TMPDIR
# ------------------------------
echo "Copying input files to TMPDIR..."

# Copy BAM + index
cp ${bam_dir}/${file} ${SLURM_TMPDIR}/
cp ${bam_dir}/${file}.bai ${SLURM_TMPDIR}/ 2>/dev/null || true

local_bam="${SLURM_TMPDIR}/${file}"

# Copy genome FASTA + index if exists
cp ${reference_dir}/GRCh38.primary_assembly.genome.fa ${SLURM_TMPDIR}/
cp ${reference_dir}/GRCh38.primary_assembly.genome.fa.fai ${SLURM_TMPDIR}/ 2>/dev/null || true

local_fasta="${SLURM_TMPDIR}/GRCh38.primary_assembly.genome.fa"

# Copy GTF
cp ${nano_dir}/liftoff_gencode/gencode.v49.annotation.gtf ${SLURM_TMPDIR}/

local_gtf="${SLURM_TMPDIR}/gencode.v49.annotation.gtf"

echo "Files in TMPDIR:"
ls -lh ${SLURM_TMPDIR}

# ------------------------------
# Run IsoQuant from TMPDIR
# ------------------------------
echo "Running IsoQuant..."

local_out="${SLURM_TMPDIR}/isoquant_${file%_restranded.sam_sorted.bam}"

python ${isoquant_dir}/isoquant.py \
    --reference ${local_fasta} \
    --data_type nanopore \
    --genedb "${local_gtf}" \
    --complete_genedb \
    -o ${local_out} \
    -t ${SLURM_CPUS_PER_TASK} \
    --sqanti_output \
    --count_exons \
    --check_canonical \
    --report_canonical all \
    --splice_correction_strategy conservative_ont \
    --matching_strategy default \
    --model_construction_strategy all \
    --report_novel_unspliced true \
    --polya_requirement never \
    --prefix MCF10a_splice_corr \
    --bam ${local_bam}

# ------------------------------
# Copy results back to SCRATCH
# ------------------------------
echo "Copying results back to SCRATCH..."
cp -r ${local_out} ${final_out_dir}/

echo "Finished at: $(date)"
