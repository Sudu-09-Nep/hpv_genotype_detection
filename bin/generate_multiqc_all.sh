#!/bin/bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

echo "═══════════════════════════════════════════════════════════"
echo "  MultiQC Report Generator"
echo "  Karkinos HPV Genotyping Pipeline"
echo "═══════════════════════════════════════════════════════════"
echo ""
echo "Project: $PROJECT_DIR"
echo ""

mkdir -p data/reports/multiqc

if ! ls data/trimmed/*_fastp.html 2>/dev/null | head -1; then
    echo "⚠️  Warning: No fastp QC reports found in data/trimmed/"
    echo "   Run the pipeline on at least one sample first."
    exit 1
fi

echo "Scanning directories for QC files:"
echo "  - data/trimmed/ (fastp reports)"
echo "  - data/alignment/ (BWA/samtools outputs)"
echo "  - data/coverage/ (coverage statistics)"
echo "  - data/logs/ (flagstat outputs)"
echo ""

# Run MultiQC (v1.35 compatible flags)
multiqc data/trimmed data/alignment data/coverage data/logs \
    --outdir data/reports/multiqc \
    --filename "hpv_pipeline_multiqc_all_samples" \
    --force

echo ""
echo "═══════════════════════════════════════════════════════════"
echo "  ✓ MultiQC Report Complete!"
echo "═══════════════════════════════════════════════════════════"
echo ""
echo "Report generated:"
echo "  HTML: data/reports/multiqc/hpv_pipeline_multiqc_all_samples.html"
echo ""
echo "To open the report:"
echo "  open data/reports/multiqc/hpv_pipeline_multiqc_all_samples.html"
echo ""
