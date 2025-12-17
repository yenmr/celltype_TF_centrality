# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(magrittr)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
n_cores <- max(25, as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE")),na.rm=T)
enableWGCNAThreads(nThreads = n_cores)

args <- commandArgs(trailingOnly = TRUE)
project      <- args[1]
rds_path_arg <- args[2]
assay        <- args[3]
celltype_col <- args[4]

# Resolve the directory of THIS script (stable even if you setwd later)
get_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", cmd, value = TRUE)
  if (length(f) == 0) return(NA_character_)
  sub("^--file=", "", f[1])
}

script_path <- get_script_path()
repo_root <- if (!is.na(script_path)) {
  # script is in <repo>/script/*.R
  normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
} else {
  # fallback: assume current working dir is repo root
  normalizePath(getwd(), mustWork = TRUE)
}
print(script_path)
print(repo_root)

# Resolve RDS path BEFORE setwd
rds_path <- normalizePath(rds_path_arg, winslash = "/", mustWork = FALSE)
if (!file.exists(rds_path)) {
  # optional convenience: if user passed only a filename, try scRNA_rds/
  candidate <- file.path(repo_root, "scRNA_rds", rds_path_arg)
  if (file.exists(candidate)) {
    rds_path <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
  } else {
    stop("Seurat RDS not found. Provided: ", rds_path_arg,
         "\nTried: ", candidate)
  }
}


dir.create(project, recursive = TRUE, showWarnings = FALSE)
setwd(project)
dir.create("kME", recursive = TRUE, showWarnings = FALSE)

seurat_obj <- readRDS(rds_path)

DefaultAssay(seurat_obj) <- assay
seurat_obj$celltype <- seurat_obj[[celltype_col]]
seurat_obj$celltype <- gsub("[^A-Za-z0-9]", "_", seurat_obj$celltype)


hsfB_list <- readr::read_tsv(file.path(repo_root, "resources", "hsfB_Osa.tsv")) %>% mutate(ID = ID_RAP)
hsfB <- hsfB_list %>% pull(ID) %>% unique()
hsfB_name <- setNames(hsfB_list$Symbol, hsfB_list$ID)
expr_mat <- LayerData(seurat_obj, assay = assay, layer = "counts")[hsfB[hsfB %in% rownames(seurat_obj)], , drop = FALSE]
n_cells_total <- ncol(seurat_obj)
hsfB_summary_overall <- tibble(
  gene          = rownames(expr_mat),
  symbol        = hsfB_name[rownames(expr_mat)],
  n_expressing  = rowSums(expr_mat > 0),
  frac_express  = n_expressing / n_cells_total
)
hsfB_summary_overall

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

# construct metacells in each group
# This step generate seurat_obj@misc$epidermal$wgcna_metacell_obj 431 metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("celltype","Orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'celltype' # set the Idents of the metacell seurat object
)
dim(seurat_obj@misc$tutorial$wgcna_metacell_obj)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)
# seurat_obj_ori <- seurat_obj
# seurat_obj <- seurat_obj_ori
celltypes <- names(table(seurat_obj$celltype))
celltype <- celltypes[1]


#saveRDS(seurat_obj, file='pre_hdWGCNA_object_root.rds')
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# This step genrate (seurat_obj@misc$epidermal$datExpr) 
# It will filter out the genes from SetupForWGCNA step and metacells from MetacellsByGroups
run_wgcna <- function(seurat_obj, celltype){
  seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = celltype, 
    group.by='celltype', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
    assay = assay, # using RNA assay
    slot = 'data' # using normalized data
  )

  # Test different soft powers:
  seurat_obj <- TestSoftPowers(
    seurat_obj,
    networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
  )
  power_table <- GetPowerTable(seurat_obj)

  # plot the results:
  plot_list <- PlotSoftPowers(seurat_obj)
  # assemble with patchwork
  p <- wrap_plots(plot_list, ncol=2)
  # ggsave(p , file = paste0(output_folder,"/SoftPower.pdf"), width = 6, height = 6)

  # construct co-expression network:
  # generate "wgcna_powerTable"   "wgcna_net"  "wgcna_modules" in names(seurat_obj@misc$epidermal)
  seurat_obj <- ConstructNetwork(
    seurat_obj,
    tom_name = celltype, # name of the topoligical overlap matrix written to disk
    overwrite_tom = TRUE
  )
  # # harmonized module eigengenes:
  # hMEs <- GetMEs(seurat_obj)
  # # module eigengenes:
  # MEs <- GetMEs(seurat_obj, harmonized=FALSE)

  seurat_obj <- ModuleEigengenes(
    seurat_obj,
    #group.by.vars="Orig.ident"
  )

  # compute eigengene-based connectivity (kME):
  seurat_obj <- ModuleConnectivity(
    seurat_obj,
    group.by = 'celltype', group_name = celltype
  )
  degree <- GetDegrees(seurat_obj)
  degree$celltype <- celltype
  degree_outfile <- file.path("kME",paste0("kME_",celltype,".tsv"))
  write_tsv(degree, degree_outfile)
  degree
  # modules_df <- GetModules(seurat_obj) %>% subset(module != 'grey')
  # table(modules_df$color)
  # modules_df$celltype <- celltype
  # modules_df
}

wgcna_list <- lapply(celltypes, function(ct) {
  message("Running WGCNA for: ", ct)

  res <- tryCatch(
    {
      df <- run_wgcna(seurat_obj, ct)
      df$celltype <- ct
      df
    },
    error = function(e) {
      message("  -> failed for ", ct, ": ", conditionMessage(e))
      NULL   # skip this cell type
    }
  )

  res
})

# wgcna_combined <- bind_rows(wgcna_list)
# write.table(wgcna_combined, file = wgcna_out_file, quote=F, row.names=F)
