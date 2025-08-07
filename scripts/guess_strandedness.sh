#!/usr/bin/env bash
# guess_standedness
# Usage: guess_strandedness.sh <genome-aligned BAM> <exon BED>

set -euo pipefail

# report usage if incorrect arguments
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <genome-aligned BAM> <exon BED>" >&2
  exit 1
fi
# assign arguments to variables
BAM=$1
BED=$2

# sanity checks
if [ ! -r "$BAM" ]; then
  echo "Error: BAM not found or unreadable: $BAM" >&2
  exit 1
fi
if [ ! -r "$BED" ]; then
  echo "Error: BED not found or unreadable: $BED" >&2
  exit 1
fi

# run RSeQC’s infer_experiment.py and capture its output
tmpf=$(mktemp)
infer_experiment.py -r "$BED" -i "$BAM" > "$tmpf"

# extract the two strandedness fractions
fwd_frac=$(grep -E 'Fraction of reads explained by "(1\+\+,1--,2\+-,2-\+|[\+][\+],[-][-])"' "$tmpf" | head -n 1 | cut -d: -f2 | xargs)
rev_frac=$(grep -E 'Fraction of reads explained by "(1\+\-,1-\+,2\+\+,2--|[\+][-],[-][\+])"' "$tmpf" | head -n 1 | cut -d: -f2 | xargs)

if [[ -z "$fwd_frac" || -z "$rev_frac" ]]; then
  echo "Error: Could not extract strandedness fractions from infer_experiment.py output." >&2
  cat "$tmpf" >&2
  rm -f "$tmpf"
  exit 1
fi

rm -f "$tmpf"

# compare and report
if awk "BEGIN { exit !($fwd_frac > $rev_frac) }"; then
  echo forward-stranded
elif awk "BEGIN { exit !($rev_frac > $fwd_frac) }"; then
  echo reverse-stranded
else
  echo unstranded
fi
