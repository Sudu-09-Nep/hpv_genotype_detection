#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params {
    project_dir = '/Users/sudarshanaryal/Desktop/test/hpv_ngs_pipeline'
    mapq        = 20
    threads     = 8
    min_length  = 36
}

process RUN_HPV_PIPELINE {

    tag "${sample_id}"

    input:
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

    R1 = '/Users/sudarshanaryal/Desktop/test/hpv_ngs_pipeline/data/raw_fastq/7936_S200_L008_R1_001.fastq'
    R2 = '/Users/sudarshanaryal/Desktop/test/hpv_ngs_pipeline/data/raw_fastq/7936_S200_L008_R2_001.fastq'

    CHUNK = tuple('7936_S200', 'PE', R1, R2)

    RUN_HPV_PIPELINE(CHUNK)

    println ""
    println "=== HPV NGS pipeline (Nextflow wrapper) finished ==="
    println "Outputs organized under data/trimmed, data/alignment, data/coverage, data/logs, data/reports"
}
