
#!/bin/bash
#SBATCH --job-name=scrna_pipeline
#SBATCH --output=pipeline_%j.out
#SBATCH --error=pipeline_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --partition=general

set -euo pipefail

# Resolve the repository root from this script's location.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd -- "${SCRIPT_DIR}/.." && pwd)}"

FASTQ_DIR="${FASTQ_DIR:-${PROJECT_DIR}/results/fastqs}"
TRANSCRIPTOME="${TRANSCRIPTOME:-${PROJECT_DIR}/references/refdata-gex-GRCh38-2024-A}"
CONDA_SH="${CONDA_SH:-}"
NXF_ENV="${NXF_ENV:-nxf_env}"
PIPELINE_MODE="${PIPELINE_MODE:-qc}"
R_SCRIPT="${R_SCRIPT:-Rscript}"
PYTHON="${PYTHON:-python}"

# Resolve relative paths before Nextflow launches tasks in work directories.
if [[ "${TRANSCRIPTOME}" != /* ]]; then
    TRANSCRIPTOME="${PROJECT_DIR}/${TRANSCRIPTOME}"
fi

# SLURM log paths are relative to the submission directory.
# Create the directory before submitting with sbatch.
cd "${PROJECT_DIR}"
mkdir -p logs

if [[ -z "${CONDA_SH}" || ! -f "${CONDA_SH}" ]]; then
    echo "ERROR: Set CONDA_SH to the actual path of conda.sh on the HPCC."
    exit 1
fi

source "${CONDA_SH}"
conda activate "${NXF_ENV}"

if ! command -v nextflow >/dev/null 2>&1; then
    echo "ERROR: Nextflow was not found in environment ${NXF_ENV}."
    exit 1
fi

if [[ ! -d "${FASTQ_DIR}" ]] ||
   ! find "${FASTQ_DIR}" -maxdepth 1 -type f -name '*.fastq.gz' -print -quit | grep -q .; then
    echo "ERROR: No .fastq.gz files found in ${FASTQ_DIR}"
    exit 1
fi

if [[ "${PIPELINE_MODE}" == "all" && ! -d "${TRANSCRIPTOME}" ]]; then
    echo "ERROR: Reference directory not found: ${TRANSCRIPTOME}"
    exit 1
fi

export CONDA_SH
export NXF_DISABLE_CHECK_LATEST=true

nextflow run "${PROJECT_DIR}/main.nf" \
    -c "${PROJECT_DIR}/nextflow.config" \
    -profile cluster \
    -resume \
    --mode "${PIPELINE_MODE}" \
    --fastq_dir "${FASTQ_DIR}" \
    --transcriptome "${TRANSCRIPTOME}" \
    --rscript "${R_SCRIPT}" \
    --python "${PYTHON}" \
    -work-dir "${PROJECT_DIR}/work_local_fastqs"
