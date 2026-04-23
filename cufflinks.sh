#!/bin/bash
#SBATCH --job-name=Cufflinks_Quantification
#SBATCH --output=cuff_all.out
#SBATCH --error=cuff_all.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --partition=medium

# 1. Environment Activation 
# activating the specialised environment containing Cufflinks v2.2.1
source /home/jvk3/.bashrc
conda activate camp_bioconda

# 2. Reference Annotation Selection 
# utilized Ensembl Release 80 (GRCh38). 
# this release is coordinate-equivalent to GENCODE v22 (Havana),
# ensuring compatibility with the Camp et al. gene models.
GTF="/scratch/alice/j/jvk3/camp_alice/ref/Homo_sapiens.GRCh38.80.gtf"
MANIFEST="/scratch/alice/j/jvk3/camp_alice/samples.txt"

# 3. Sequential Processing Loop 
# loop iterates through the full 730-sample manifest.
for FULL_PATH in $(cat $MANIFEST)
do
    SRR_ID=$(basename $FULL_PATH)
    
    echo "Processing: $SRR_ID ---"
    
    # path definitions for input (TopHat BAM) and output (Cufflinks Results)
    BAM="/scratch/alice/j/jvk3/camp_alice/results/${SRR_ID}_tophat/accepted_hits.bam"
    OUT_DIR="/scratch/alice/j/jvk3/camp_alice/results/${SRR_ID}_cufflinks_V2"
    
    # check for BAM existence before executing
    if [ -f "$BAM" ]; then
        mkdir -p $OUT_DIR
        
        # running cufflinks:
        # -p 8: multi-threaded processing
        # -G: quantify based on supplied GTF only (quantification mode, not assembly)
        cufflinks -p 8 -G $GTF -o $OUT_DIR $BAM
        
        echo "Successfully quantified $SRR_ID"
    else
        echo "WARNING: Skipping $SRR_ID. BAM file not found."
    fi
done

