#!/usr/bin/env bash
set -euo pipefail

# Destination base directory
BASE=/ceph/cms/store/user/tvami/DiMuonPlusX/Ntuples_v1.0.1/

# ------------------------------------------------------------------
# 1) Create full directory tree
# ------------------------------------------------------------------
for sample in Signal Data BkgMC; do
  for region in sr vr; do
    mkdir -p "$BASE/$sample/$region"
  done
done

# ------------------------------------------------------------------
# 2) Move skimmed ROOT files from ./ into the tree
# ------------------------------------------------------------------
shopt -s nullglob
for f in ./skimmed_*.root ./trigger_study_*.root; do
  fname=$(basename "$f")

  # Determine output base: append "_trigger_study" for trigger_study files
  if [[ "$fname" == trigger_study_* ]]; then
    OUTBASE="${BASE%/}_trigger_study/"
    prefix=trigger_study
  else
    OUTBASE="$BASE"
    prefix=skimmed
  fi

  # Extract region
  if [[ "$fname" == ${prefix}_sr_* ]]; then
    region=sr
  elif [[ "$fname" == ${prefix}_vr_* ]]; then
    region=vr
  else
    echo "WARNING: cannot determine region for $fname — skipping"
    continue
  fi

  # Determine sample type
  if [[ "$fname" == *TpTp* ]]; then
    sample=Signal
  elif [[ "$fname" == *ScoutingPFRun3* ]]; then
    sample=Data
  else
    sample=BkgMC
  fi

  mkdir -p "$OUTBASE/$sample/$region"
  mv "$f" "$OUTBASE/$sample/$region/"
done
