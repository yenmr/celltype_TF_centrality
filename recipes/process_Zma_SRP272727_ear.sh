#!/usr/bin/bash
Rscript --slave script/run_hdWGCNA_Zma.R Zma_SRP272727_ear SRP272727_23_26.rds SCT Celltype
Rscript --slave script/combine_kME_file.R Zma_SRP272727_ear
Rscript --slave script/analyze_coexpression_Zma.R Zma_SRP272727_ear

