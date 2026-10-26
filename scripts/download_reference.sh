
#!/usr/bin/env bash
#SBATCH --job-name=download_ref
#SBATCH --output=download_ref_%j.out
#SBATCH --error=download_ref_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=general

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_DIR}/references}"
REFERENCE_NAME="refdata-gex-GRCh38-2024-A"
ARCHIVE="${REFERENCE_DIR}/${REFERENCE_NAME}.tar.gz"
REFERENCE_PATH="${REFERENCE_DIR}/${REFERENCE_NAME}"
URL="https://cf.10xgenomics.com/supp/cell-exp/${REFERENCE_NAME}.tar.gz"

mkdir -p "${REFERENCE_DIR}"

if [[ -d "${REFERENCE_PATH}" && -s "${REFERENCE_PATH}/reference.json" ]]; then
    echo "Reference already exists and reference.json is present:"
    echo "${REFERENCE_PATH}"
    exit 0
fi

if [[ ! -s "${ARCHIVE}" ]]; then
    echo "Downloading 10x Genomics reference: ${REFERENCE_NAME}"
    curl --fail --location --retry 3 --retry-delay 5 \
        --output "${ARCHIVE}.part" "${URL}"
    mv "${ARCHIVE}.part" "${ARCHIVE}"
else
    echo "Using existing archive: ${ARCHIVE}"
fi

if ! gzip -t "${ARCHIVE}"; then
    echo "ERROR: The reference archive failed gzip validation: ${ARCHIVE}" >&2
    exit 1
fi

echo "Validating archive contents..."
if ! tar -tzf "${ARCHIVE}" | grep -Fxq "${REFERENCE_NAME}/reference.json"; then
    echo "ERROR: Archive does not contain the expected ${REFERENCE_NAME}/reference.json." >&2
    exit 1
fi

echo "Extracting reference..."
tar -xzf "${ARCHIVE}" -C "${REFERENCE_DIR}"

if [[ ! -s "${REFERENCE_PATH}/reference.json" ]]; then
    echo "ERROR: Extraction did not produce the expected reference.json." >&2
    exit 1
fi

echo "Reference download and extraction completed: ${REFERENCE_PATH}"
