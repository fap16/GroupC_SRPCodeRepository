#!/usr/bin/bash
#SBATCH --job-name=sra_to_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-16%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aakg1@student.le.ac.uk
#SBATCH --output=/scratch/alice/a/aakg1/SRA_project/logs/sra_to_fastq.out
#SBATCH --error=/scratch/alice/a/aakg1/SRA_project/logs/sra_to_fastq.err
#SBATCH --export=NONE

# Load SRA toolkit
module load sratoolkit/3.0.0-5fetwpi

# Assigning directories for use
PROJ="/scratch/alice/a/aakg1/SRA_project"
ACC_LIST="${PROJ}/SRR_Acc_List.txt"
SRA_DIR="${PROJ}/sra"
FASTQ_DIR="${PROJ}/fastq"
LOG_DIR="${PROJ}/logs"

# check if accession list is present, if present the code will run
[[ -f "${ACC_LIST}" ]] || { echo "ERROR: Accession list not found: ${ACC_LIST}"; exit 1; }

cd "${PROJ}"

# Creating required directories
mkdir -p "${SRA_DIR}" "${FASTQ_DIR}" "${LOG_DIR}"

# loop for task array to run through all accessions
for i in $(seq "${SLURM_ARRAY_TASK_ID}" "${WORKERS}" "${TOTAL}"); do
    srr=$(sed -n "${i}p" "${ACC_LIST}" | tr -d '[:space:]')
    [[ -z "${srr}" ]] && continue

    # creating a temporary directory to avoid overwriting or corrupting files
    TMP_DIR="${PROJ}/tmp/${SLURM_ARRAY_JOB_ID}_${i}"
    mkdir -p "${TMP_DIR}"

    # Downloading all srr
    prefetch --output-directory "${SRA_DIR}" "${srr}" \

  { echo "ERROR: prefetch failed for ${srr}"; continue; } #will show which srr failed 

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

done
