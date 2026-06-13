#!/bin/bash
set -euo pipefail

SAMPLE="$1"
MODE="$2"
R1="$3"
R2="${4:-}"
REFERENCE="/Users/sudarshanaryal/Desktop/test/hpv_ngs_pipeline/data/reference/reference.fasta"
MAPQ_THRESHOLD="${5:-20}"
THREADS="${6:-8}"
MIN_LENGTH="${7:-36}"

if [ "$MODE" = "SE" ]; then
    R2=""
fi

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

mkdir -p data/trimmed data/alignment data/coverage data/reports data/logs

echo "═══════════════════════════════════════════════════════════"
echo "  Karkinos HPV Genotyping Pipeline"
echo "  Sample: $SAMPLE"
echo "  Mode: $MODE"
echo "  Reference: $(basename "$REFERENCE")"
echo "  MAPQ Threshold: $MAPQ_THRESHOLD"
echo "  Threads: $THREADS"
echo "  Min Read Length: $MIN_LENGTH"
echo "═══════════════════════════════════════════════════════════"

echo "[1/6] Running fastp QC..."
if [ "$MODE" = "PE" ]; then
  fastp --thread "$THREADS" \
    --detect_adapter_for_pe \
    --overrepresentation_analysis \
    --correction \
    --cut_right \
    --length_required "$MIN_LENGTH" \
    --html "data/trimmed/${SAMPLE}_fastp.html" \
    --json "data/trimmed/${SAMPLE}_fastp.json" \
    --report_title "HPV Sample $SAMPLE QC Report" \
    -i "$R1" -I "$R2" \
    -o "data/trimmed/${SAMPLE}_R1_paired.fastq" \
    -O "data/trimmed/${SAMPLE}_R2_paired.fastq" \
    --unpaired1 "data/trimmed/${SAMPLE}_R1_unpaired.fastq" \
    --unpaired2 "data/trimmed/${SAMPLE}_R2_unpaired.fastq"
else
  fastp --thread "$THREADS" \
    --overrepresentation_analysis \
    --correction \
    --cut_right \
    --length_required "$MIN_LENGTH" \
    --html "data/trimmed/${SAMPLE}_fastp.html" \
    --json "data/trimmed/${SAMPLE}_fastp.json" \
    --report_title "HPV Sample $SAMPLE QC Report" \
    -i "$R1" \
    -o "data/trimmed/${SAMPLE}_paired.fastq"
fi
echo "✓ fastp complete: data/trimmed/${SAMPLE}_fastp.html"

echo "[2/6] Preparing reference genome..."
if [[ ! -f "${REFERENCE}.bwt" ]]; then
  echo "  Indexing reference genome with bwa..."
  bwa index "$REFERENCE"
fi
echo "✓ Reference ready"

echo "[3/6] Aligning to reference..."
if [ "$MODE" = "PE" ]; then
  bwa mem -t "$THREADS" "$REFERENCE" \
    "data/trimmed/${SAMPLE}_R1_paired.fastq" \
    "data/trimmed/${SAMPLE}_R2_paired.fastq" \
    | samtools sort -@ "$THREADS" -o "data/alignment/${SAMPLE}.sorted.bam"
else
  bwa mem -t "$THREADS" "$REFERENCE" \
    "data/trimmed/${SAMPLE}_paired.fastq" \
    | samtools sort -@ "$THREADS" -o "data/alignment/${SAMPLE}.sorted.bam"
fi
samtools index "data/alignment/${SAMPLE}.sorted.bam"
echo "✓ Alignment complete: data/alignment/${SAMPLE}.sorted.bam"

echo "[4/6] Filtering by MAPQ ≥ $MAPQ_THRESHOLD..."
samtools view -b -q "$MAPQ_THRESHOLD" "data/alignment/${SAMPLE}.sorted.bam" \
  > "data/alignment/${SAMPLE}.filtered_MAPQ${MAPQ_THRESHOLD}.bam"
samtools index "data/alignment/${SAMPLE}.filtered_MAPQ${MAPQ_THRESHOLD}.bam"
echo "✓ Filtering complete"

echo "[5/6] Generating alignment statistics..."
samtools idxstats "data/alignment/${SAMPLE}.filtered_MAPQ${MAPQ_THRESHOLD}.bam" \
  > "data/alignment/${SAMPLE}_idxstats.txt"
samtools flagstat "data/alignment/${SAMPLE}.filtered_MAPQ${MAPQ_THRESHOLD}.bam" \
  > "data/logs/${SAMPLE}_flagstat.txt"
samtools coverage "data/alignment/${SAMPLE}.filtered_MAPQ${MAPQ_THRESHOLD}.bam" \
  > "data/coverage/${SAMPLE}_coverage.txt"
echo "✓ Statistics complete"

echo "[6/6] Generating HPV genotype report..."

REPORT_FILE="data/reports/${SAMPLE}_genotype_report.txt"
FLAGSTAT_FILE="data/logs/${SAMPLE}_flagstat.txt"
COVERAGE_FILE="data/coverage/${SAMPLE}_coverage.txt"
IDXSTATS_FILE="data/alignment/${SAMPLE}_idxstats.txt"

set +e

ALIGNED="NA"
UNMAPPED="NA"
if [ -f "$FLAGSTAT_FILE" ]; then
  ALIGNED=$(grep "reads mapped" "$FLAGSTAT_FILE" 2>/dev/null | head -1 | awk '{print $1}')
  UNMAPPED=$(grep "reads unmapped" "$FLAGSTAT_FILE" 2>/dev/null | head -1 | awk '{print $1}')
fi

{
  echo "═══════════════════════════════════════════════════════════"
  echo "  Karkinos HPV Genotyping Report"
  echo "═══════════════════════════════════════════════════════════"
  echo ""
  echo "Sample: $SAMPLE"
  echo "Mode: $MODE"
  echo "Reference: $(basename "$REFERENCE")"
  echo "MAPQ Threshold: $MAPQ_THRESHOLD"
  echo "Threads: $THREADS"
  echo "Min Read Length: $MIN_LENGTH"
  echo ""
  echo "───────────────────────────────────────────────────────────"
  echo "  Alignment Statistics"
  echo "───────────────────────────────────────────────────────────"
  echo ""

  if [ -f "$FLAGSTAT_FILE" ]; then
    echo "Flagstat summary:"
    cat "$FLAGSTAT_FILE"
  else
    echo "Flagstat summary: (file not found: $FLAGSTAT_FILE)"
  fi
  echo ""

  if [ -f "$COVERAGE_FILE" ]; then
    echo "Coverage summary:"
    cat "$COVERAGE_FILE"
  else
    echo "Coverage summary: (file not found: $COVERAGE_FILE)"
  fi
  echo ""

  if [ -f "$IDXSTATS_FILE" ]; then
    echo "Idxstats (per-chromosome mapping):"
    cat "$IDXSTATS_FILE"
  else
    echo "Idxstats: (file not found: $IDXSTATS_FILE)"
  fi
  echo ""
  echo "───────────────────────────────────────────────────────────"
  echo "  Genotype Call"
  echo "───────────────────────────────────────────────────────────"
  echo ""

  REF_NAME=$(basename "$REFERENCE" | sed 's/\.fasta$//' | sed 's/\.fa$//')
  echo "Detected genotype: $REF_NAME (based on reference)"
  echo ""
  echo "Total aligned reads: ${ALIGNED:-NA}"
  echo "Unmapped reads: ${UNMAPPED:-NA}"
  echo ""
  echo "═══════════════════════════════════════════════════════════"
  echo "  End of HPV Genotyping Report"
  echo "═══════════════════════════════════════════════════════════"
} > "$REPORT_FILE"

echo "✓ Genotype report: $REPORT_FILE"

set -e
