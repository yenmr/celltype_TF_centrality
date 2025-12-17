library(dplyr)
library(readr)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
project <- args[1] # "Mpo_snRNA7d_tak1"
setwd(project)

source("../lib/lib_ortho_GO.R")
go_gene_ontology_table <- getQueryGOFromAth("Ath")

cytokinin <- read_tsv(gzfile("../resources/Cytokinin-golden-list.tsv.gz"))

hsfB_list <- read_tsv("../resources/hsfB_Ath.tsv")# %>% mutate(ID = ID_RAP)
hsfB <- hsfB_list %>% pull(ID) %>% unique()
hsf_name <- setNames(hsfB_list$Symbol, hsfB_list$ID)

run_analysis <- function(df, target_hsf) {
  print(paste("Now working on", target_hsf))
  df_name <- deparse(substitute(df))
  hsf_module <- df %>% filter(gene_name %in% target_hsf) %>%
    rename(kME = weighted_degree) %>%
    dplyr::select(celltype, module, kME, tf_rank, tf_total, tf_score)

  module_total_count <- df %>%
    count(module, celltype) %>%
    semi_join(hsf_module, by = c("celltype", "module"))

  module_cytokinin_count <- df %>% 
    filter(gene_name %in% cytokinin$ID) %>%
    count(module, celltype) %>%
    semi_join(hsf_module, by = c("celltype", "module"))

  total_count <- df %>% pull(gene_name) %>% unique() %>% length()

  cytokinin_count <- df %>% filter(gene_name %in% cytokinin$ID) %>%
    pull(gene_name) %>% unique() %>% length()

  df_stats <- hsf_module %>% filter(module != "grey") %>%
    inner_join(module_cytokinin_count,
      by = c("module", "celltype")
    ) %>% rename(n_cyt = n) %>%
    inner_join(
      module_total_count,
      by = c("module", "celltype")
    ) %>% rename(n_total = n) %>%
    rowwise() %>%
    mutate(
      a = n_cyt,
      b = n_total - n_cyt,
      c = cytokinin_count - a,
      d = total_count - a - b - c,
      tbl = list(matrix(c(a, b, c, d), nrow = 2)),
      ft  = list(fisher.test(tbl, alternative = "greater")),
      OR  = ft$estimate,
      pvalue   = ft$p.value
    ) %>%
    ungroup() %>%
    mutate(qvalue = p.adjust(pvalue, method = "BH")) %>%
    arrange(pvalue) %>%                     # sort by q
    rename(
      CK_mod     = a,
      CK_out_mod      = b,
      nonCK_mod  = c,
      nonCK_out_mod   = d,
      rank_kME_TF   = tf_rank,
      nTF_module = tf_total,
      rank_pct_TF   = tf_score,
      p_value    = pvalue,
      q_value    = qvalue,
      kME_query  = kME
    ) %>% mutate(query_tf = hsf_name[target_hsf]) %>%
    dplyr::select(celltype, module, CK_mod, CK_out_mod, nonCK_mod, nonCK_out_mod, OR, p_value, q_value, query_tf, kME_query, rank_kME_TF, nTF_module, rank_pct_TF)
    
  OR_cytokinin_outfile <- paste0(prefix,"_",hsf_name[target_hsf],"_OR.tsv")
  write_tsv(df_stats, OR_cytokinin_outfile)

  #GO analysis
  bg <- df$gene_name %>% unique()
  df_GO <- df %>% semi_join(hsf_module, by = c("celltype", "module"))

  GO_results <- df_GO %>%
    group_by(celltype) %>%
    group_split() %>%
    purrr::map_df(function(df_sub) {
      print(df_sub$celltype %>% unique())
      genes <- df_sub$gene_name
      tb <- run_GO(genes, bg)
      if (is.null(tb)) {return(NULL)}
      tb %>% mutate(celltype = unique(df_sub$celltype))
    })
  out_GO_file <- paste0(prefix,"_",hsf_name[target_hsf],"_GO.tsv")
  write_tsv(GO_results, out_GO_file)
}

root <- read_tsv("kME_combined.tsv")
prefix <- project

tf <- read_tsv(gzfile("../resources/Ath_TF_iTAK_list.txt.gz"), col_names=FALSE)  %>%
  mutate(Gene_ID = sub("\\.\\d+$","",X1), TF_family = X2) %>% select(Gene_ID, TF_family) %>% distinct(Gene_ID, .keep_all=TRUE)
tf_info <- root %>% 
  inner_join(tf, by = c("gene_name" = "Gene_ID")) %>%
  group_by(module, celltype) %>%
  arrange(desc(degree), .by_group = TRUE) %>%
  mutate(
    tf_rank = row_number(),
    tf_total = n(),
    tf_score = ifelse(tf_total > 1, (tf_total - tf_rank) / (tf_total - 1), 1)
  ) %>%
  ungroup() 
write_tsv(tf_info,paste0(prefix, "_TF_rank.tsv"))
root <- root %>% left_join(tf_info %>% select(-degree, -weighted_degree), by = c("gene_name", "module", "celltype")) 

hsfB_list %>% filter(ID %in% root$gene_name) %>% pull(ID) %>% purrr::set_names() %>% purrr::map(~ run_analysis(root, .x))
# run_analysis(root, "AT4G36990")
