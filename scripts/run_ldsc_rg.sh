#!/usr/bin/env bash
set -euo pipefail

: "${LDSC_PY:?Set LDSC_PY to the path of ldsc.py}"
: "${LDSC_LD_PREFIX:?Set LDSC_LD_PREFIX to the reference LD-score prefix}"
: "${LDSC_WEIGHT_PREFIX:?Set LDSC_WEIGHT_PREFIX to the weight prefix}"
: "${HM3_SNPLIST:?Set HM3_SNPLIST to w_hm3.snplist}"

left=${1:?Usage: run_ldsc_rg.sh LEFT_SUMSTATS RIGHT_SUMSTATS OUTPUT_PREFIX}
right=${2:?Usage: run_ldsc_rg.sh LEFT_SUMSTATS RIGHT_SUMSTATS OUTPUT_PREFIX}
output=${3:?Usage: run_ldsc_rg.sh LEFT_SUMSTATS RIGHT_SUMSTATS OUTPUT_PREFIX}

python3 "$LDSC_PY" \
  --rg "$left,$right" \
  --ref-ld-chr "$LDSC_LD_PREFIX" \
  --w-ld-chr "$LDSC_WEIGHT_PREFIX" \
  --print-snps "$HM3_SNPLIST" \
  --print-cov \
  --out "$output"
