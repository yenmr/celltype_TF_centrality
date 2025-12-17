#!/usr/bin/bash
Rscript --slave script/run_hdWGCNA_Zma.R Zma_SRP335180_root SRP335180.rds SCT Celltype
Rscript --slave script/combine_kME_file.R Zma_SRP335180_root
Rscript --slave script/analyze_coexpression_Zma.R Zma_SRP335180_root

