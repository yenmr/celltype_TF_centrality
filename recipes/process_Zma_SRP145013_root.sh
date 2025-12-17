#!/usr/bin/bash
Rscript --slave script/run_hdWGCNA_Zma.R Zma_SRP145013_root SRP145013.rds SCT Celltype
Rscript --slave script/combine_kME_file.R Zma_SRP145013_root
Rscript --slave script/analyze_coexpression_Zma.R Zma_SRP145013_root

