#!/bin/bash
#SBATCH --job-name=Camp_Loop
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/alice/j/jvk3/camp_alice_scratch/logs/loop.log

module load sratoolkit
OUT_DIR="/scratch/alice/j/jvk3/camp_alice/fastq"
mkdir -p /scratch/alice/j/jvk3/sra_cache
export VDB_CONFIG_DIR=/scratch/alice/j/jvk3/sra_cache
vdb-config --set /repository/user/main/public/root=/scratch/alice/j/jvk3/sra_cache

# Loop reads list starting from line 65
sed -n '65,734p' ~/camp_alice/meta/srr_list.txt | while read SRR; do
    if [ ! -f "${OUT_DIR}/${SRR}_1.fastq" ]; then
        echo "Starting $SRR..."
        fasterq-dump "$SRR" --outdir "$OUT_DIR" --threads 4 --split-files
        rm -rf /scratch/alice/j/jvk3/sra_cache/sra/${SRR}.sra
    fi
done
