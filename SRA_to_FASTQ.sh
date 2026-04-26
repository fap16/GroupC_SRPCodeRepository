#Slurm script for SRR accession conversion to FASTQ
#last update 27/04/2026

#!/usr/bin/bash
#SBATCH --job-name=sra_to_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-16%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user= (input personal email)
#SBATCH --output=/scratch/alice/a/aakg1/SRA_project/logs/sra_to_fastq_%A_%a.out
#SBATCH --error=/scratch/alice/a/aakg1/SRA_project/logs/sra_to_fastq_%A_%a.err
#SBATCH --export=NONE

# Load SRA toolkit
module load sratoolkit/3.0.0-5fetwpi

# Assigning directories for use
PROJ="/scratch/alice/a/aakg1/SRA_project"
ACC_LIST="${PROJ}/SRR_Acc_List.txt"
SRA_DIR="${PROJ}/sra"
FASTQ_DIR="${PROJ}/fastq"
LOG_DIR="${PROJ}/logs"
TMP_all="${PROJ}/tmp"

# set total job size and workers needed for loop
TOTAL=734
WORKERS=16

# check if accession list is present, if present the code will run
[[ -f "${ACC_LIST}" ]] || { echo "ERROR: Accession list not found: ${ACC_LIST}"; exit 1; }
cd "${PROJ}"

# Creating required directories
mkdir -p "${SRA_DIR}" "${FASTQ_DIR}" "${LOG_DIR}" "${TMP_all}"

# loop for task array to run through all accessions
for i in $(seq "${SLURM_ARRAY_TASK_ID}" "${WORKERS}" "${TOTAL}"); do
    srr=$(sed -n "${i}p" "${ACC_LIST}" | tr -d '[:space:]')
    [[ -z "${srr}" ]] && continue

#store files in temporary to avoid overwriting
TMP_DIR="${TMP_all}/${SLURM_ARRAY_JOB_ID}_${i}"
mkdir -p "${TMP_DIR}"

 # Prefetching/Downloading all srr
prefetch --output-directory "${SRA_DIR}" "${srr}" \
|| { echo "ERROR: prefetch failed for ${srr}"; continue; } #will show which srr failed 
SRA_PATH="${SRA_DIR}/${srr}/${srr}.sra"
    [[ -f "${SRA_PATH}" ]] || { echo "ERROR: SRA file missing for ${srr}"; continue; }

# Convert to FASTQ for STAR alignment
fasterq-dump \
        --split-files \
        --threads "${SLURM_CPUS_PER_TASK:-8}" \
        --outdir "${FASTQ_DIR}" \
        --temp "${TMP_DIR}" \
        "${SRA_PATH}" \
        || { echo "ERROR: fasterq-dump failed for ${srr}"; continue; }

# Compress FASTQ files
 gzip -f "${FASTQ_DIR}/${srr}"*.fastq

rm -rf "${TMP_DIR}" #delete temporary file folder

    echo "Success: ${srr} completed successfully"
done
