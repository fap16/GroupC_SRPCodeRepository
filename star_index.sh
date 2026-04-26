#!/usr/bin/bash
#SBATCH --job-name=STAR_index
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --mail-user=aakg1@student.le.ac.uk
#SBATCH --mail-type=END,BEGIN,FAIL
#SBATCH --output=/scratch/alice/a/aakg1/genome_files/star_index_%j.out
#SBATCH --error=/scratch/alice/a/aakg1/genome_files/star_index_%j.err

set -euo pipefail

module load star/2.7.11a-5kstrrw
#set own paths for use and directories generated

GENOME_DIR="/scratch/alice/a/aakg1/genome_files/star_index"
FASTA="/scratch/alice/a/aakg1/genome_files/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF="/scratch/alice/a/aakg1/genome_files/Homo_sapiens.GRCh38.110.gtf"

mkdir -p "$GENOME_DIR"

STAR \
--runThreadN ${SLURM_CPUS_PER_TASK} \
--runMode genomeGenerate \
--genomeDir "$GENOME_DIR" \
--genomeFastaFiles "$FASTA" \
--sjdbGTFfile "$GTF" \
--sjdbOverhang 100
