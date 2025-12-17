#!/usr/bin/bash
Rscript --slave script/run_hdWGCNA_Zma.R Zma_SRP224648_leaf SRP224648.rds SCT Celltype
Rscript --slave script/combine_kME_file.R Zma_SRP224648_leaf
Rscript --slave script/analyze_coexpression_Zma.R Zma_SRP224648_leaf

