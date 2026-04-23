#!/bin/bash
#SBATCH --job-name=TopHat_Replication_Batch
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00

# 1. Environment & Legacy Compatibility 
# TopHat 2.1.2 runs with Python 2.7. 
# call the binary directly via absolute path to ensure version consistency.
module load bowtie2/2.5.1
source activate camp_bioconda
PY2="/home/j/jvk3/.conda/envs/camp_bioconda/bin/python2.7"
TOPHAT_BIN="/cm/shared/spack/opt/spack/linux-rocky9-x86_64_v3/gcc-12.3.0/tophat-2.1.2-ay5hc7g6ermzbdxb3uolul7dbeaurixh/bin/tophat"

# 2. Reference & Annotation
# using GENCODE v22 (No-Chr) to match Ensembl-style GRCh38 indexing
MANIFEST="/scratch/alice/j/jvk3/camp_alice/samples.txt"
INDEX="/scratch/alice/j/jvk3/camp_alice/ref/GRCh38_ensembl"
GTF="/scratch/alice/j/jvk3/camp_alice/ref/gencode.v22.no_chr.gtf"

# 3. Parallel Batch Processing 
# processing a specific slice of the 730 samples (in this batch 1-15)
# allows for manageable job sizes on the cluster.
for i in {1..15}
do
    # use 'sed' to pull the specific sample path from manifest
    SAMPLE_BASE=$(sed -n "${i}p" $MANIFEST)
    if [ -z "$SAMPLE_BASE" ]; then continue; fi
    CELL=$(basename $SAMPLE_BASE)
    
    echo "Processing Sample $i: $CELL"

    # Alignment Execution:
    # --no-novel-juncs ensures strict replication of the gene models in the paper.
    $PY2 $TOPHAT_BIN -p 8 \
        -G $GTF \
        --no-novel-juncs \
        -o /scratch/alice/j/jvk3/camp_alice/results/${CELL}_tophat \
        $INDEX \
        ${SAMPLE_BASE}_1.fastq ${SAMPLE_BASE}_2.fastq
done
