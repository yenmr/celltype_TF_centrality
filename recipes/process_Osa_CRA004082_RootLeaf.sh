#!/usr/bin/bash
Rscript --slave script/run_hdWGCNA_Osa.R Osa_CRA004082_RootLeaf CRA004082.rds RNA Celltype
Rscript --slave script/combine_kME_file.R Osa_CRA004082_RootLeaf
Rscript --slave script/analyze_coexpression_Osa.R Osa_CRA004082_RootLeaf

