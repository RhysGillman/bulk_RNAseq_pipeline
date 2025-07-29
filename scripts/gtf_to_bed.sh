#!/bin/bash
#gtf_to_bed

in="$1"
out="$2"

awk '
  $3=="exon" {
    # extract gene_id from the 9th column
    match($0, /gene_id "([^"]+)"/, m)
    gene=m[1]
    # GTF is 1-based inclusive; BED uses 0-based [start,end)
    printf("%s\t%d\t%d\t%s\t.\t%s\n",
      $1, $4-1, $5, gene, $7)
  }
' "$in" \
  | sort -k1,1 -k2,2n \
  > "$out"
