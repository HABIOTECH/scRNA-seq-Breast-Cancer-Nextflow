
#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
OUTDIR="${OUTDIR:-${PROJECT_DIR}/results/fastqs}"
MAX_PARALLEL="${MAX_PARALLEL:-2}"

mkdir -p "${OUTDIR}" "${PROJECT_DIR}/logs"

download_one() {
    local sample="$1"
    local srr="$2"
    local lane="$3"
    local work="${OUTDIR}/${srr}"

    local out1="${OUTDIR}/${sample}_S1_L${lane}_R1_001.fastq.gz"
    local out2="${OUTDIR}/${sample}_S1_L${lane}_R2_001.fastq.gz"

    if [[ -s "${out1}" && -s "${out2}" ]] &&
       gzip -t "${out1}" "${out2}" 2>/dev/null; then
        echo "SKIP: ${srr} already complete"
        return 0
    fi

    mkdir -p "${work}"

    local report
    if ! report="$(wget -qO- \
        "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srr}&result=read_run&fields=run_accession,fastq_ftp&format=tsv")"; then
        echo "ERROR: Could not retrieve ENA report for ${srr}" >&2
        return 1
    fi

    local -a urls
    mapfile -t urls < <(
        printf '%s\n' "${report}" |
        awk -F '\t' 'NR==2 {print $2}' |
        tr ';' '\n' |
        sed '/^[[:space:]]*$/d'
    )

    if [[ "${#urls[@]}" -ne 2 ]]; then
        echo "ERROR: Expected two FASTQ URLs for ${srr}" >&2
        return 1
    fi

    local url
    for url in "${urls[@]}"; do
        if ! wget -c --tries=20 --waitretry=30 --timeout=120 \
            -P "${work}" "https://${url}"; then
            echo "ERROR: Download failed for ${srr}: ${url}" >&2
            return 1
        fi
    done

    local r1="${work}/${srr}_1.fastq.gz"
    local r2="${work}/${srr}_2.fastq.gz"

    if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
        echo "ERROR: Missing or empty FASTQ file for ${srr}" >&2
        return 1
    fi

    if ! gzip -t "${r1}" "${r2}"; then
        echo "ERROR: FASTQ gzip validation failed for ${srr}" >&2
        return 1
    fi

    local r1_len r2_len
    r1_len="$(zcat "${r1}" | awk 'NR==2 {print length; exit}')"
    r2_len="$(zcat "${r2}" | awk 'NR==2 {print length; exit}')"

    if [[ "${r1_len}" != "28" || ! "${r2_len}" =~ ^[0-9]+$ ]] ||
       (( r2_len <= 30 )); then
        echo "ERROR: Unexpected read lengths for ${srr}: R1=${r1_len:-unknown} R2=${r2_len:-unknown}" >&2
        return 1
    fi

    mv "${r1}" "${out1}"
    mv "${r2}" "${out2}"
    rmdir "${work}" 2>/dev/null || true

    echo "DONE: ${srr} -> ${sample}_S1_L${lane}"
}

RUNS=(
"hMG3 SRR31298425 001"
"hMG2 SRR31298426 001"
"hMG1 SRR31298427 001"
"BBM5 SRR31298428 001"
"BBM5 SRR31298429 002"
"BBM5 SRR31298430 003"
"BBM5 SRR31298431 004"
"BBM5 SRR31298432 005"
"BSM1 SRR31298433 001"
"BSM1 SRR31298434 002"
"BSM1 SRR31298435 003"
"BSM1 SRR31298436 004"
"BSM1 SRR31298437 005"
"BSM1 SRR31298438 006"
"BSM1 SRR31298439 007"
"BSM1 SRR31298440 008"
"BBM0 SRR31298441 001"
"BBM0 SRR31298442 002"
"BBM0 SRR31298443 003"
"BBM0 SRR31298444 004"
"BBM0 SRR31298445 005"
"BBM0 SRR31298446 006"
"BBM0 SRR31298447 007"
"BBM0 SRR31298448 008"
)

pids=()
failed=0

for row in "${RUNS[@]}"; do
    read -r sample srr lane <<< "${row}"

    download_one "${sample}" "${srr}" "${lane}" \
        > "${PROJECT_DIR}/logs/${srr}.log" 2>&1 &
    pids+=("$!")

    if [[ "${#pids[@]}" -ge "${MAX_PARALLEL}" ]]; then
        for pid in "${pids[@]}"; do
            if ! wait "${pid}"; then
                echo "ERROR: A download task failed (PID ${pid})." >&2
                failed=1
            fi
        done
        pids=()
    fi
done

for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
        echo "ERROR: A download task failed (PID ${pid})." >&2
        failed=1
    fi
done

if (( failed != 0 )); then
    echo "ERROR: One or more downloads failed. Check logs in ${PROJECT_DIR}/logs." >&2
    exit 1
fi

echo "SUCCESS: All FASTQ downloads completed or were already validated."
