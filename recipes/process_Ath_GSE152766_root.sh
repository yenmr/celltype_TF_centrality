#!/usr/bin/bash
Rscript --slave script/run_hdWGCNA_Ath.R Ath_GSE152766_root GSE152766_Root_Atlas.rds RNA time.celltype.anno
Rscript --slave script/combine_kME_file.R Ath_GSE152766_root
Rscript --slave script/analyze_coexpression_Ath.R Ath_GSE152766_root

