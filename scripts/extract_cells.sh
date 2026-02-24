#!/bin/bash
#SBATCH --job-name=Camp_Extract
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/extract_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jvk3@le.ac.uk

module load sra-toolkit/3.0.0

ACC_LIST="meta/srr_list.txt"
OUT_DIR="fastq"

mkdir -p $OUT_DIR

while read -r SRR; do
    echo "Processing $SRR"
    prefetch "$SRR"
    # Tell the tool to use more threads since now have 48 requested
    fasterq-dump "$SRR" --outdir "$OUT_DIR" --threads 12 --split-files
done < "$ACC_LIST"
