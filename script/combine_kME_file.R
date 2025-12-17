#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

# pattern for your files
project <- args[1] # "Mpo_snRNA7d_tak1"
setwd(project)

files <- list.files(path = "kME", pattern = "^kME_.*\\.tsv$",full.names = TRUE)
print(files)
message("Found ", length(files), " files")

# read and add celltype column
kme_list <- lapply(files, function(f) {
  message("Reading: ", f)
  
  df <- read_tsv(f, show_col_types = FALSE)
  
  # extract celltype from filename: kME_root_XXX.tsv -> XXX
  celltype <- basename(f) %>%
    str_remove("^kME_root_") %>%
    str_remove("\\.tsv$")
  
  df %>%
    mutate(celltype = celltype)
})

# combine all
kme_combined <- bind_rows(kme_list)

# optional: inspect
print(dim(kme_combined))
print(head(kme_combined))

# write to file
out_file <- "kME_combined.tsv"
write_tsv(kme_combined, out_file)
message("Combined file written to: ", out_file)

