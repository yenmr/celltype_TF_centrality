#!/usr/bin/bash
Rscript --slave script/run_hdWGCNA_Ath.R Ath_GSE297576_seedling GSE297576_seurat_object.thaliana_seedling_atlas.RDS SCT celltype
Rscript --slave script/combine_kME_file.R Ath_GSE297576_seedling
Rscript --slave script/analyze_coexpression_Ath.R Ath_GSE297576_seedling

