#!/usr/bin/bash
Rscript --slave script/run_hdWGCNA_Sly.R Sly_SRP280242_shoot SRP280242.rds RNA Celltype
Rscript --slave script/combine_kME_file.R Sly_SRP280242_shoot
Rscript --slave script/analyze_coexpression_Sly.R Sly_SRP280242_shoot

