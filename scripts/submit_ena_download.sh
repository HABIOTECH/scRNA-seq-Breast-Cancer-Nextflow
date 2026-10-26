
#!/usr/bin/env bash
#SBATCH --job-name=ena_fastq
#SBATCH --output=logs/ena_%j.out
#SBATCH --error=logs/ena_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=168:00:00
#SBATCH --partition=general

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd -- "${SCRIPT_DIR}/.." && pwd)}"

cd "${PROJECT_DIR}"
mkdir -p logs

OUTDIR="${OUTDIR:-${PROJECT_DIR}/results/fastqs}"
export OUTDIR

echo "Project directory: ${PROJECT_DIR}"
echo "FASTQ output directory: ${OUTDIR}"

bash "${SCRIPT_DIR}/download_ena_fastq.sh"
