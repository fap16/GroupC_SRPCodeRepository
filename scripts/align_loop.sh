#!/bin/bash
#SBATCH --job-name=Align_Loop
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/alice/j/jvk3/camp_alice_scratch/logs/align_loop.log

module load bowtie2
module load samtools

FASTQ_DIR="/scratch/alice/j/jvk3/camp_alice/fastq"
INDEX="/home/j/jvk3/camp_alice/ref/hg19_index"
OUT_DIR="/scratch/alice/j/jvk3/camp_alice/aligned"
mkdir -p $OUT_DIR

# This loop looks at your list and only processes SRRs that have FASTQ files ready
sed -n '1,734p' ~/camp_alice/meta/srr_list.txt | while read SRR; do
    # Only run if the FASTQ exists AND the sorted BAM doesn't exist yet
    if [ -f "${FASTQ_DIR}/${SRR}_1.fastq" ] && [ ! -f "${OUT_DIR}/${SRR}_sorted.bam" ]; then
        echo "Processing $SRR..."
        
        # Align
        bowtie2 -p 8 -q -x $INDEX -1 ${FASTQ_DIR}/${SRR}_1.fastq -2 ${FASTQ_DIR}/${SRR}_2.fastq -S ${OUT_DIR}/${SRR}.sam
        
        # Convert, Sort, and Clean
        samtools view -bS ${OUT_DIR}/${SRR}.sam > ${OUT_DIR}/${SRR}.bam
        samtools sort ${OUT_DIR}/${SRR}.bam -o ${OUT_DIR}/${SRR}_sorted.bam
        rm ${OUT_DIR}/${SRR}.sam ${OUT_DIR}/${SRR}.bam
        
        echo "Finished $SRR."
    fi
done
