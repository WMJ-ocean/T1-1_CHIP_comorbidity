#!/bin/bash
# 23_run_sens.sh —— 敏感性 MR 批处理
cd "$(dirname "$0")"
R="/d/R-4.5.2/bin/x64/Rscript.exe"
SENS="../05_iv/sens"
OUTDIR="../07_mr/sens"
mkdir -p "$OUTDIR"

log() { echo "[$(date +%H:%M:%S)] $*"; }

# ---------- A. 正向：CHIP 胚系 各敏感集 → 五病 ----------
for SET in noTERT noMHC; do
  for P in MPN PV MDS AA ITP; do
    IV="$SENS/CHIP_germline_${SET}.tsv"
    TGT="../06_targets/${P}_targets.tsv"
    TAG="CHIPg${SET}_to_${P}"
    log "RUN $TAG"
    "$R" 20_mr.R "$IV" "$TGT" "$TAG" > "$OUTDIR/${TAG}.log" 2>&1
    log "  done $TAG (exit=$?)"
  done
done

# ---------- B. 反向：MPN/PV 去 JAK2 / 去 JAK2+TERT → CHIP ----------
for SET in noJAK2 noJAK2TERT; do
  for P in MPN PV; do
    IV="$SENS/${P}_${SET}.tsv"
    TGT="../06_targets/CHIP_targets.tsv"
    TAG="${P}${SET}_to_CHIP"
    log "RUN $TAG"
    "$R" 20_mr.R "$IV" "$TGT" "$TAG" > "$OUTDIR/${TAG}.log" 2>&1
    log "  done $TAG (exit=$?)"
  done
done

log "ALL DONE"
