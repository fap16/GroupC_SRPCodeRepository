#!/usr/bin/bash
#SBATCH --job-name=sra_to_fastq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-16%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --output=/scratch/alice/j/jvk3/camp_alice/logs/sra_to_fastq.out
#SBATCH --error=/scratch/alice/j/jvk3/camp_alice/logs/sra_to_fastq.err
#SBATCH --export=NONE

# Load SRA toolkit
module load sratoolkit/3.0.0-5fetwpi

# Assigning directories for use/can set your own

PROJ="/scratch/alice/j/jvk3/camp_alice"
ACC_LIST="${PROJ}/SRR_Acc_List.txt"
SRA_DIR="${PROJ}/sra"
FASTQ_DIR="${PROJ}/fastq"
LOG_DIR="${PROJ}/logs"
TMP_BASE="${PROJ}/tmp"

TOTAL=734 #set your own total or how many you need
WORKERS=16 #how many jobs running/limit is 16 for AlICE

# Check accession list exists
[[ -f "${ACC_LIST}" ]] || { echo "ERROR: Accession list not found: ${ACC_LIST}"; exit 1; }
cd "${PROJ}" 

# Create required directories
mkdir -p "${SRA_DIR}" "${FASTQ_DIR}" "${LOG_DIR}" "${TMP_BASE}"

# Loop through accessions assigned to this array task
for i in $(seq "${SLURM_ARRAY_TASK_ID}" "${WORKERS}" "${TOTAL}"); do
    srr=$(sed -n "${i}p" "${ACC_LIST}" | tr -d '[:space:]')
    [[ -z "${srr}" ]] && continue
    
    # Temporary directory for this accession
    TMP_DIR="${TMP_BASE}/${SLURM_ARRAY_JOB_ID}_${i}"
    mkdir -p "${TMP_DIR}"
   
    # Download SRA
prefetch --output-directory "${SRA_DIR}" "${srr}" \
|| { echo "ERROR: prefetch failed for ${srr}"; rm -rf "${TMP_DIR}"; continue; }

# Check SRA exists
SRA_PATH="${SRA_DIR}/${srr}/${srr}.sra"
    [[ -f "${SRA_PATH}" ]] || { echo "ERROR: SRA file missing for ${srr}"; continue; }
    
    # Convert to FASTQ
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
