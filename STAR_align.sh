#STAR ALIGNMENT code
# last update 26/04/2026

#!/bin/bash
#SBATCH --job-name=STAR_alignment
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/scratch/alice/a/aakg1/SRA_project/logs/star_%j.out
#SBATCH --error=/scratch/alice/a/aakg1/SRA_project/logs/star_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aakg1@student.le.ac.uk

# Load modules
module load star/2.7.11a-5kstrrw
module load samtools/1.17-wenuvv5

# paths/use own paths
GENOME_DIR="/scratch/alice/a/aakg1/genome_files/star_index" 
FASTQ_DIR="/scratch/alice/a/aakg1/SRA_project/fastq"
OUT_DIR="/scratch/alice/a/aakg1/SRA_project/alignment_results"

#generate output directory

mkdir -p $OUT_DIR

# Loop  all R1 files
for R1 in ${FASTQ_DIR}/*_1.fastq.gz; do
# identify 2nd pair
R2=${R1/_1.fastq.gz/_2.fastq.gz}
# base naming
SAMPLE=$(basename $R1 _1.fastq.gz)
 # 4. Run STAR alignment using genome index
    STAR --runThreadN 16 \
         --genomeDir $GENOME_DIR \
         --readFilesIn $R1 $R2 \
         --readFilesCommand zcat \
         --outFileNamePrefix ${OUT_DIR}/${SAMPLE}_ \
         --outSAMtype BAM SortedByCoordinate \ #setting the parameters
         --quantMode GeneCounts \
         --outSAMunmapped Within \
         --limitBAMsortRAM 30000000000 #indicates how much memory allowed to use in RAM(bytes)

 # 5. Index the BAM
    samtools index ${OUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam

#end
done

echo "All samples complete!"
