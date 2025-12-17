#!/usr/bin/bash
Rscript --slave script/run_hdWGCNA_Mpo.R Mpo_snRNA7d_tak1 tak1_snRNA.rds RNA labels
Rscript --slave script/combine_kME_file.R Mpo_snRNA7d_tak1
Rscript --slave script/analyze_coexpression_Mpo.R Mpo_snRNA7d_tak1

