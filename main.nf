nextflow.enable.dsl=2

/*
 * Run modes:
 *   --mode qc    Use the existing FASTQs, run FastQC, and create a raw-read
 *                MultiQC report. Stop here for review.
 *   --mode all   Run the complete workflow, including Cell Ranger and the
 *                final MultiQC report.
 *
 * Input FASTQs are already downloaded. Files are selected from params.fastq_dir
 * using the naming convention:
 *   SAMPLE_S1_L###_R1_001.fastq.gz
 *   SAMPLE_S1_L###_R2_001.fastq.gz
 *
 * Sequencing lanes are deliberately NOT concatenated. Files from all lanes
 * belonging to one biological sample are grouped and passed together to
 * Cell Ranger. The L### names are retained so Cell Ranger can identify and
 * combine the lanes correctly.
 *
 * These are 10x Chromium 3' v3.1 libraries. R1 contains the 28-bp barcode/UMI
 * read, so trimming is intentionally not performed.
 */


process FASTQC {
    tag "${sample}"
    publishDir "${params.outdir}/fastqc", mode: 'copy'
    cpus 4
    memory '8 GB'
    time '4h'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("*_fastqc.zip"), emit: qc_zips

    script:
    """
    set -euo pipefail
    module load ${params.fastqc_module}
    fastqc --threads ${task.cpus} --noextract ${reads}
    """
}

process MULTIQC_RAW {
    tag 'raw-read-qc'
    conda "${projectDir}/envs/multiqc_env.yml"
    publishDir "${params.outdir}/multiqc/raw", mode: 'copy'
    cpus 2
    memory '4 GB'
    time '2h'

    input:
    path qc_zips

    output:
    path 'multiqc_raw.html', emit: report

    script:
    """
    set -euo pipefail
    multiqc --version | grep -F "1.35" >/dev/null || {
        echo "ERROR: Expected MultiQC 1.35."
        multiqc --version
        exit 1
    }
    multiqc --force --filename multiqc_raw.html .
    """
}

process CELLRANGER_COUNT {
    tag "${sample}"
    publishDir "${params.outdir}/cellranger", mode: 'copy'
    cpus 16
    memory '128 GB'
    time '48h'

    input:
    tuple val(sample), path(reads), path(raw_qc_report)

    output:
    tuple val(sample), path("${sample}/outs/filtered_feature_bc_matrix"), emit: count_matrix
    tuple val(sample), path("${sample}/metrics_summary.csv"), emit: metrics
    tuple val(sample), path("${sample}/web_summary.html"), emit: web_summary

    script:
    """
    set -euo pipefail
    module load cellranger/10.0.0
    mkdir -p tmp
    export TMPDIR="\$PWD/tmp"
    test -s ${raw_qc_report}

    # Keep the original lane-specific filenames. Symlinks avoid copying the
    # large input FASTQs into the Cell Ranger work directory.
    mkdir -p input_fastqs
    for fq in ${reads}; do
        ln -s "\$(readlink -f "\$fq")" "input_fastqs/\$(basename "\$fq")"
    done

    if [[ ! -d "${params.transcriptome}" ]]; then
        echo "ERROR: Cell Ranger reference directory not found: ${params.transcriptome}"
        exit 1
    fi

    cellranger count --id=${sample} \
                     --transcriptome=${params.transcriptome} \
                     --fastqs=input_fastqs \
                     --sample=${sample} \
                     --create-bam=false \
                     --localcores=${task.cpus} \
                     --localmem=${task.memory.toGiga()}

    cp ${sample}/outs/metrics_summary.csv ${sample}/metrics_summary.csv
    cp ${sample}/outs/web_summary.html ${sample}/web_summary.html
    """
}

process MULTIQC_FINAL {
    tag 'final-qc'
    conda "${projectDir}/envs/multiqc_env.yml"
    publishDir "${params.outdir}/multiqc/final", mode: 'copy'
    cpus 2
    memory '8 GB'
    time '2h'

    input:
    path fastqc_zips
    path cellranger_metrics, stageAs: 'cellranger_metrics_*'

    output:
    path 'multiqc_final.html', emit: report

    script:
    """
    set -euo pipefail
    multiqc --version | grep -F "1.35" >/dev/null || {
        echo "ERROR: Expected MultiQC 1.35."
        multiqc --version
        exit 1
    }
    multiqc . --force --filename multiqc_final.html
    """
}

process SEURAT_ANALYSIS {
    tag "${sample}"
    conda "${projectDir}/envs/seurat_env.yml"
    publishDir "${params.outdir}/seurat", mode: 'copy'
    cpus 4
    memory '16 GB'
    time '4h'

    input:
    tuple val(sample), path(matrix)

    output:
    tuple val(sample), path("${sample}_seurat.pdf"), emit: plots
    tuple val(sample), path("${sample}_seurat.rds"), emit: objects

    script:
    """
    set -euo pipefail
    Rscript "${projectDir}/bin/seurat_script.R" --input "${matrix}" --output "${sample}_seurat"
    """
}

process SCANPY_ANALYSIS {
    tag "${sample}"
    conda "${projectDir}/envs/scanpy_env.yml"
    publishDir "${params.outdir}/scanpy", mode: 'copy'
    cpus 4
    memory '16 GB'
    time '4h'

    input:
    tuple val(sample), path(matrix)

    output:
    tuple val(sample), path("${sample}_scanpy.pdf"), emit: plots
    tuple val(sample), path("${sample}_scanpy.h5ad"), emit: objects

    script:
    """
    set -euo pipefail
    python "${projectDir}/bin/scanpy_script.py" --input "${matrix}" --output "${sample}_scanpy"
    """
}

workflow {
    if (!(params.mode in ['qc', 'all'])) {
        error "Invalid --mode '${params.mode}'. Use --mode qc or --mode all."
    }

    // Parse every downloaded FASTQ and group all lanes by biological sample.
    // The regex also validates that files follow the expected 10x naming scheme.
    grouped_fastqs_ch = Channel
        .fromPath("${params.fastq_dir}/*.fastq.gz", checkIfExists: true)
        .map { file ->
            def match = file.name =~ /^(.+)_S1_(L[0123456789][0123456789][0123456789])_R([12])_001[.]fastq[.]gz$/
            if (!match.matches()) {
                throw new IllegalArgumentException(
                    "Unexpected FASTQ filename: ${file.name}. Expected SAMPLE_S1_L###_R[12]_001.fastq.gz"
                )
            }
            tuple(match.group(1), file)
        }
        .groupTuple()
        .map { sample, files ->
            def sorted_files = files.sort { a, b -> a.name <=> b.name }
            def r1_lanes = sorted_files
                .findAll { it.name.contains('_R1_') }
                .collect { (it.name =~ /_S1_(L[0123456789][0123456789][0123456789])_R1_001[.]fastq[.]gz$/)[0][1] }
                .toSet()
            def r2_lanes = sorted_files
                .findAll { it.name.contains('_R2_') }
                .collect { (it.name =~ /_S1_(L[0123456789][0123456789][0123456789])_R2_001[.]fastq[.]gz$/)[0][1] }
                .toSet()

            if (r1_lanes != r2_lanes) {
                throw new IllegalArgumentException(
                    "Unpaired lanes for sample ${sample}: R1=${r1_lanes}, R2=${r2_lanes}"
                )
            }
            tuple(sample, sorted_files)
        }

    FASTQC(grouped_fastqs_ch)

    fastqc_files_ch = FASTQC.out.qc_zips
        .map { sample, files -> files }
        .flatten()
        .collect()

    MULTIQC_RAW(fastqc_files_ch)

    if (params.mode == 'all') {
        // The raw MultiQC report is an explicit dependency of Cell Ranger.
        mapping_ch = grouped_fastqs_ch
            .combine(MULTIQC_RAW.out.report)
            .map { sample, reads, raw_report ->
                tuple(sample, reads, raw_report)
            }

        CELLRANGER_COUNT(mapping_ch)

        SEURAT_ANALYSIS(CELLRANGER_COUNT.out.count_matrix)
        SCANPY_ANALYSIS(CELLRANGER_COUNT.out.count_matrix)

        metrics_files_ch = CELLRANGER_COUNT.out.metrics
            .map { sample, file -> file }
            .flatten()
            .collect()

        MULTIQC_FINAL(fastqc_files_ch, metrics_files_ch)
    }
}
