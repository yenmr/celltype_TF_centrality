#!/usr/bin/bash
Rscript --slave script/run_hdWGCNA_Sly.R Sly_SRP286427_root SRP286427.rds SCT Celltype
Rscript --slave script/combine_kME_file.R Sly_SRP286427_root
Rscript --slave script/analyze_coexpression_Sly.R Sly_SRP286427_root

