#!/bin/bash
set -euo pipefail

# Usage:
#   PE:
#     ./pipeline.sh SAMPLE PE R1.fastq R2.fastq [MAPQ] [THREADS] [MIN_LENGTH] [MIN_READS] [MIN_ABUNDANCE_PCT] [MIN_MEAN_DEPTH] [AMPLICON_BED] [MIN_AMPLICON_COVERAGE_PCT]
#   SE:
#     ./pipeline.sh SAMPLE SE R1.fastq "" [MAPQ] [THREADS] [MIN_LENGTH] [MIN_READS] [MIN_ABUNDANCE_PCT] [MIN_MEAN_DEPTH] [AMPLICON_BED] [MIN_AMPLICON_COVERAGE_PCT]
#
# Example:
#   ./pipeline.sh S1 PE R1.fastq R2.fastq 30 8 36 100 0.1 20 data/reference/amplicons.bed 80
#
# Notes:
#   - AMPLICON_BED is optional. If provided, contig names in the BED must match the FASTA headers/reference names.
#   - If AMPLICON_BED is not provided, the amplicon coverage filter is reported as NA and is not used to fail calls.

SAMPLE="${1:?ERROR: SAMPLE is required}"
MODE="${2:?ERROR: MODE is required: PE or SE}"
R1="${3:?ERROR: R1 FASTQ is required}"
R2="${4:-}"

REFERENCE="/Users/sudarshanaryal/Desktop/test/hpv_ngs_pipeline/data/reference/reference.fasta"
MAPQ_THRESHOLD="${5:-20}"
THREADS="${6:-8}"
MIN_LENGTH="${7:-36}"

# Genotype-calling thresholds
MIN_READS="${8:-100}"
MIN_ABUNDANCE_PCT="${9:-0.1}"
MIN_MEAN_DEPTH="${10:-20}"
AMPLICON_BED="${11:-}"
MIN_AMPLICON_COVERAGE_PCT="${12:-80}"

if [ "$MODE" = "SE" ]; then
  R2=""
elif [ "$MODE" != "PE" ]; then
  echo "ERROR: MODE must be PE or SE" >&2
  exit 1
fi

if [ "$MODE" = "PE" ] && [ -z "$R2" ]; then
  echo "ERROR: R2 FASTQ is required for PE mode" >&2
  exit 1
fi

if [ -n "$AMPLICON_BED" ] && [ ! -f "$AMPLICON_BED" ]; then
  echo "ERROR: AMPLICON_BED not found: $AMPLICON_BED" >&2
  exit 1
fi

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

mkdir -p data/trimmed data/alignment data/coverage data/reports data/logs data/calls

FILTERED_BAM="data/alignment/${SAMPLE}.filtered_MAPQ${MAPQ_THRESHOLD}.bam"
FLAGSTAT_FILE="data/logs/${SAMPLE}_flagstat.txt"
COVERAGE_FILE="data/coverage/${SAMPLE}_coverage.txt"
IDXSTATS_FILE="data/alignment/${SAMPLE}_idxstats.txt"
AMPLICON_DEPTH_FILE="data/coverage/${SAMPLE}_amplicon_depth.txt"
AMPLICON_SUMMARY_FILE="data/coverage/${SAMPLE}_amplicon_summary.tsv"
CALLS_FILE="data/calls/${SAMPLE}_genotype_calls.tsv"
CALLS_TMP_FILE="data/calls/${SAMPLE}_genotype_calls.unsorted.tsv"
REPORT_FILE="data/reports/${SAMPLE}_genotype_report.txt"

MIN_MEAN_MAPQ="$MAPQ_THRESHOLD"

printf '%s\n' "═══════════════════════════════════════════════════════════"
printf '%s\n' "  Karkinos HPV Genotyping Pipeline"
printf '%s\n' "  Sample: $SAMPLE"
printf '%s\n' "  Mode: $MODE"
printf '%s\n' "  Reference: $(basename "$REFERENCE")"
printf '%s\n' "  MAPQ Threshold: $MAPQ_THRESHOLD"
printf '%s\n' "  Threads: $THREADS"
printf '%s\n' "  Min Read Length: $MIN_LENGTH"
printf '%s\n' "  Min Genotype Reads: $MIN_READS"
printf '%s\n' "  Min Relative Abundance: ${MIN_ABUNDANCE_PCT}%"
printf '%s\n' "  Min Mean Depth: ${MIN_MEAN_DEPTH}x"
printf '%s\n' "  Min Mean MAPQ: $MIN_MEAN_MAPQ"
if [ -n "$AMPLICON_BED" ]; then
  printf '%s\n' "  Amplicon BED: $AMPLICON_BED"
  printf '%s\n' "  Min Amplicon Coverage: ${MIN_AMPLICON_COVERAGE_PCT}%"
else
  printf '%s\n' "  Amplicon BED: not provided; amplicon coverage filter will be NA"
fi
printf '%s\n' "═══════════════════════════════════════════════════════════"

echo "[1/7] Running fastp QC..."
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

echo "[2/7] Preparing reference genome..."
if [[ ! -f "${REFERENCE}.bwt" ]]; then
  echo "  Indexing reference genome with bwa..."
  bwa index "$REFERENCE"
fi
if [[ ! -f "${REFERENCE}.fai" ]]; then
  echo "  Indexing reference genome with samtools faidx..."
  samtools faidx "$REFERENCE"
fi
echo "✓ Reference ready"

echo "[3/7] Aligning to reference..."
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

echo "[4/7] Filtering by MAPQ ≥ $MAPQ_THRESHOLD..."
samtools view -b -q "$MAPQ_THRESHOLD" "data/alignment/${SAMPLE}.sorted.bam" > "$FILTERED_BAM"
samtools index "$FILTERED_BAM"
echo "✓ Filtering complete: $FILTERED_BAM"

echo "[5/7] Generating alignment statistics..."
samtools idxstats "$FILTERED_BAM" > "$IDXSTATS_FILE"
samtools flagstat "$FILTERED_BAM" > "$FLAGSTAT_FILE"
samtools coverage "$FILTERED_BAM" > "$COVERAGE_FILE"
echo "✓ Statistics complete"

echo "[6/7] Calculating amplicon coverage and genotype calls..."

if [ -n "$AMPLICON_BED" ]; then
  # Per-base depth across the amplicon BED. BED coordinates are 0-based half-open;
  # samtools depth outputs 1-based positions.
  samtools depth -a -b "$AMPLICON_BED" "$FILTERED_BAM" > "$AMPLICON_DEPTH_FILE"

  awk 'BEGIN {
         OFS="\t";
         print "rname", "amplicon_bases", "covered_bases", "amplicon_coverage_pct", "amplicon_mean_depth";
       }
       {
         r=$1; d=$3 + 0;
         total[r]++;
         depthsum[r]+=d;
         if (d > 0) covered[r]++;
       }
       END {
         for (r in total) {
           cov=(total[r] > 0 ? covered[r] / total[r] * 100 : 0);
           md=(total[r] > 0 ? depthsum[r] / total[r] : 0);
           printf "%s\t%d\t%d\t%.4f\t%.4f\n", r, total[r], covered[r], cov, md;
         }
       }' "$AMPLICON_DEPTH_FILE" | sort -k1,1 > "$AMPLICON_SUMMARY_FILE"
else
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "rname" "amplicon_bases" "covered_bases" "amplicon_coverage_pct" "amplicon_mean_depth" \
    > "$AMPLICON_SUMMARY_FILE"
fi

# Build final genotype-call table from samtools coverage.
# samtools coverage columns:
#   rname startpos endpos numreads covbases coverage meandepth meanbaseq meanmapq
awk -v min_reads="$MIN_READS" \
    -v min_abund="$MIN_ABUNDANCE_PCT" \
    -v min_depth="$MIN_MEAN_DEPTH" \
    -v min_mapq="$MIN_MEAN_MAPQ" \
    -v amp_file="$AMPLICON_SUMMARY_FILE" \
    -v amp_bed="$AMPLICON_BED" \
    -v min_amp_cov="$MIN_AMPLICON_COVERAGE_PCT" \
    'BEGIN {
       FS=OFS="\t";
       has_amp=(amp_bed != "");
       while ((getline line < amp_file) > 0) {
         if (line ~ /^rname\t/) continue;
         split(line, a, "\t");
         amp_cov[a[1]]=a[4] + 0;
         amp_md[a[1]]=a[5] + 0;
       }
       close(amp_file);
     }
     /^#/ { next }
     NF >= 9 {
       rname=$1;
       numreads[rname]=$4 + 0;
       ref_coverage[rname]=$6 + 0;
       meandepth[rname]=$7 + 0;
       meanbaseq[rname]=$8 + 0;
       meanmapq[rname]=$9 + 0;
       total += numreads[rname];
       order[++n]=rname;
     }
     END {
       print "genotype", "numreads", "abundance_pct", "reference_coverage_pct", "meandepth", "meanbaseq", "meanmapq", "amplicon_coverage_pct", "amplicon_mean_depth", "PASS_READS", "PASS_ABUNDANCE", "PASS_MEAN_DEPTH", "PASS_MEAN_MAPQ", "PASS_AMPLICON_COVERAGE", "FINAL_CALL";

       for (i=1; i<=n; i++) {
         r=order[i];
         abund=(total > 0 ? numreads[r] / total * 100 : 0);

         pass_reads=(numreads[r] >= min_reads ? "PASS" : "FAIL");
         pass_abund=(abund >= min_abund ? "PASS" : "FAIL");
         pass_depth=(meandepth[r] >= min_depth ? "PASS" : "FAIL");
         pass_mapq=(meanmapq[r] >= min_mapq ? "PASS" : "FAIL");

         if (has_amp) {
           acov=(r in amp_cov ? amp_cov[r] : 0);
           amd=(r in amp_md ? amp_md[r] : 0);
           pass_amp=(acov >= min_amp_cov ? "PASS" : "FAIL");
         } else {
           acov="NA";
           amd="NA";
           pass_amp="NA";
         }

         final="FAIL";
         if (pass_reads == "PASS" && pass_abund == "PASS" && pass_depth == "PASS" && pass_mapq == "PASS" && (!has_amp || pass_amp == "PASS")) {
           final="PASS";
         }

         if (has_amp) {
           printf "%s\t%d\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\t%s\t%s\t%s\t%s\n", r, numreads[r], abund, ref_coverage[r], meandepth[r], meanbaseq[r], meanmapq[r], acov, amd, pass_reads, pass_abund, pass_depth, pass_mapq, pass_amp, final;
         } else {
           printf "%s\t%d\t%.6f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", r, numreads[r], abund, ref_coverage[r], meandepth[r], meanbaseq[r], meanmapq[r], acov, amd, pass_reads, pass_abund, pass_depth, pass_mapq, pass_amp, final;
         }
       }
     }' "$COVERAGE_FILE" > "$CALLS_TMP_FILE"

{
  head -n 1 "$CALLS_TMP_FILE"
  tail -n +2 "$CALLS_TMP_FILE" | sort -t $'\t' -k15,15r -k3,3nr
} > "$CALLS_FILE"
rm -f "$CALLS_TMP_FILE"

echo "✓ Genotype calls: $CALLS_FILE"
if [ -n "$AMPLICON_BED" ]; then
  echo "✓ Amplicon coverage: $AMPLICON_SUMMARY_FILE"
fi

echo "[7/7] Generating HPV genotype report..."

set +e

ALIGNED="NA"
UNMAPPED="NA"
if [ -f "$FLAGSTAT_FILE" ]; then
  ALIGNED=$(grep "reads mapped" "$FLAGSTAT_FILE" 2>/dev/null | head -1 | awk '{print $1}')
  UNMAPPED=$(grep "reads unmapped" "$FLAGSTAT_FILE" 2>/dev/null | head -1 | awk '{print $1}')
fi

PASSING_GENOTYPES=$(awk -F '\t' 'NR > 1 && $15 == "PASS" {print $1}' "$CALLS_FILE" | paste -sd "," -)
if [ -z "$PASSING_GENOTYPES" ]; then
  PASSING_GENOTYPES="None"
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
  echo "Genotype-calling thresholds:"
  echo "  Minimum genotype reads: $MIN_READS"
  echo "  Minimum relative abundance: ${MIN_ABUNDANCE_PCT}%"
  echo "  Minimum mean depth: ${MIN_MEAN_DEPTH}x"
  echo "  Minimum mean MAPQ: $MIN_MEAN_MAPQ"
  if [ -n "$AMPLICON_BED" ]; then
    echo "  Amplicon BED: $AMPLICON_BED"
    echo "  Minimum amplicon coverage: ${MIN_AMPLICON_COVERAGE_PCT}%"
  else
    echo "  Amplicon BED: not provided"
    echo "  Amplicon coverage filter: NA / not applied"
  fi
  echo ""
  echo "───────────────────────────────────────────────────────────"
  echo "  Alignment Statistics"
  echo "───────────────────────────────────────────────────────────"
  echo ""

  if [ -f "$FLAGSTAT_FILE" ]; then
    echo "Flagstat summary:"
    cat "$FLAGSTAT_FILE"
  else
    echo "Flagstat summary: file not found: $FLAGSTAT_FILE"
  fi
  echo ""

  echo "Total aligned reads: ${ALIGNED:-NA}"
  echo "Unmapped reads: ${UNMAPPED:-NA}"
  echo ""

  if [ -f "$COVERAGE_FILE" ]; then
    echo "Raw samtools coverage summary:"
    cat "$COVERAGE_FILE"
  else
    echo "Coverage summary: file not found: $COVERAGE_FILE"
  fi
  echo ""

  if [ -f "$IDXSTATS_FILE" ]; then
    echo "Idxstats per-reference mapping:"
    cat "$IDXSTATS_FILE"
  else
    echo "Idxstats: file not found: $IDXSTATS_FILE"
  fi
  echo ""

  if [ -n "$AMPLICON_BED" ] && [ -f "$AMPLICON_SUMMARY_FILE" ]; then
    echo "Amplicon coverage summary:"
    cat "$AMPLICON_SUMMARY_FILE"
    echo ""
  fi

  echo "───────────────────────────────────────────────────────────"
  echo "  Genotype Calls"
  echo "───────────────────────────────────────────────────────────"
  echo ""
  echo "Passing genotype(s): $PASSING_GENOTYPES"
  echo ""
  echo "Filtered genotype table:"
  if [ -f "$CALLS_FILE" ]; then
    cat "$CALLS_FILE"
  else
    echo "Genotype-call table not found: $CALLS_FILE"
  fi
  echo ""
  echo "Interpretation:"
  echo "  A genotype is called only if it passes read-count, abundance, mean-depth, and mean-MAPQ filters."
  if [ -n "$AMPLICON_BED" ]; then
    echo "  Amplicon coverage is also required because an amplicon BED was provided."
  else
    echo "  Amplicon coverage is not required because no amplicon BED was provided."
  fi
  echo ""
  echo "═══════════════════════════════════════════════════════════"
  echo "  End of HPV Genotyping Report"
  echo "═══════════════════════════════════════════════════════════"
} > "$REPORT_FILE"

echo "✓ Genotype report: $REPORT_FILE"

set -e

echo "Done."
