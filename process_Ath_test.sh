#!/usr/bin/bash
Rscript --slave script/run_hdWGCNA_Ath.R Ath_test Ath_test.rds SCT celltype
Rscript --slave script/combine_kME_file.R Ath_test
Rscript --slave script/analyze_coexpression_Ath.R Ath_test

#bash run_pipeline.sh --species Ath --dataset Ath_test --rds scRNA_rds/Ath_test.rds --assay SCT --celltype celltype
