#!/usr/bin/env bash
set -euo pipefail

# usage:
# bash run_pipeline.sh --species Ath --dataset Ath_test --rds scRNA_rds/Ath_test.rds --assay SCT --celltype celltype

while [[ $# -gt 0 ]]; do
  case "$1" in
    --species)  SPECIES="$2"; shift 2 ;;
    --dataset)  DATASET_DIR="$2"; shift 2 ;;
    --rds)      SEURAT_RDS="$2"; shift 2 ;;
    --assay)    ASSAY="$2"; shift 2 ;;
    --celltype) CELLTYPE_COL="$2"; shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

: "${SPECIES:?Missing --species}"
: "${DATASET_DIR:?Missing --dataset}"
: "${SEURAT_RDS:?Missing --rds}"
: "${ASSAY:?Missing --assay}"
: "${CELLTYPE_COL:?Missing --celltype}"

case "$SPECIES" in
  Ath|Mpo|Osa|Sly|Zma) ;;
  *) echo "Unsupported species: $SPECIES (use Ath/Mpo/Osa/Sly/Zma)" >&2; exit 1 ;;
esac

RUN_SCRIPT="script/run_hdWGCNA_${SPECIES}.R"
ANA_SCRIPT="script/analyze_coexpression_${SPECIES}.R"

echo "[RUN] Rscript --slave $RUN_SCRIPT $DATASET_DIR $SEURAT_RDS $ASSAY $CELLTYPE_COL"
Rscript --slave "$RUN_SCRIPT" "$DATASET_DIR" "$SEURAT_RDS" "$ASSAY" "$CELLTYPE_COL"

echo "[RUN] Rscript --slave script/combine_kME_file.R $DATASET_DIR"
Rscript --slave script/combine_kME_file.R "$DATASET_DIR"

echo "[RUN] Rscript --slave $ANA_SCRIPT $DATASET_DIR"
Rscript --slave "$ANA_SCRIPT" "$DATASET_DIR"

