#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params {
    project_dir  = '/Users/sudarshanaryal/Desktop/test/hpv_ngs_pipeline'
    mapq         = 20
    threads      = 8
    min_length   = 36
    samplesheet  = '${project_dir}/data/samplesheet.csv'
}

process RUN_HPV_PIPELINE {

    tag "${sample_id}"

    input:
    // Use val() so Nextflow treats these as plain strings (full paths), not staged files
    tuple val(sample_id), val(mode), val(r1_str), val(r2_str)

    script:
    """
    echo "=== Running Karkinos HPV pipeline for sample: ${sample_id} ==="

    cd "${params.project_dir}"

    chmod +x bin/run_hpv_pipeline.sh

    ./bin/run_hpv_pipeline.sh \
        "${sample_id}" \
        "${mode}" \
        "${r1_str}" \
        "${r2_str}" \
        "${params.mapq}" \
        "${params.threads}" \
        "${params.min_length}"
    """
}

workflow {

    samples_ch = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)
        .map { row ->
            tuple(
                row.sample_id,
                row.mode,
                row.fastq_r1,  // full path from CSV
                row.fastq_r2   // full path from CSV
            )
        }

    RUN_HPV_PIPELINE(samples_ch)

    println ""
    println "=== HPV NGS pipeline (Nextflow wrapper) finished ==="
    println "Outputs organized under data/trimmed, data/alignment, data/coverage, data/logs, data/reports"
}
