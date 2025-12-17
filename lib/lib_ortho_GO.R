lib_file <- sys.frame(1)$ofile
lib_dir <- dirname(normalizePath(lib_file))
ortho_file <- file.path(lib_dir, "Orthogroups_long.tsv.gz")


get_go_term2gene <- function(ont = c("BP", "MF", "CC")) {
  ont <- match.arg(ont)

  go_data <- clusterProfiler:::get_GO_data(
    "org.At.tair.db",  # OrgDb
    ont,               # ontology
    "TAIR"             # keytype
  )
  gs <- DOSE:::getGeneSet(go_data)
  tibble::tibble(
    GO       = rep(names(gs), lengths(gs)),
    TAIR     = unlist(gs, use.names = FALSE),
    ONTOLOGY = ont
  )
}

protein2gene_rules <- list(
  default = function(x) sub("\\..*$", "", x),
  Ath = function(x) sub("\\..*$", "", x),
  Mpo = function(x) sub("\\..*$", "", x),
  Sly = function(x) sub("T(\\d+)\\.\\d+", "G\\1", x),
  Osa = function(x) sub("\\..*$", "", x),
  Zma = function(x) sub("_P\\d+$", "", x)
)

gene_from_protein <- function(protein, species) {
  f <- protein2gene_rules[[species]]
  if (is.null(f)) f <- protein2gene_rules$default
  f(protein)
}

getQueryGOFromAth <- function(query_species) {
  go_annotations_subject <- dplyr::bind_rows(
    get_go_term2gene("BP"),
    get_go_term2gene("MF"),
    get_go_term2gene("CC")
  )

  go_terms <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = unique(na.omit(go_annotations_subject$GO)),
    keytype = "GOID",
    columns = "TERM"
  ) |> tibble::as_tibble()

  go_annotations_subject <- go_annotations_subject |> 
    dplyr::left_join(go_terms, by = c("GO" = "GOID"))

  if (query_species == "Ath") {
	  return(go_annotations_subject |> dplyr::mutate(Gene = TAIR) |> dplyr::select(GO, Gene, ONTOLOGY, TERM))
  }

  ortho <- read.delim(gzfile(ortho_file)) |>
    dplyr::filter(Species %in% c("Ath", query_species)) |>
    tibble::as_tibble() |> dplyr::mutate(Gene = purrr::map2_chr(Protein, Species, ~ gene_from_protein(.x, .y)))

  ortho_Ath <- ortho |> dplyr::filter(Species =='Ath')
  ortho_query <- ortho |> dplyr::filter(Species == query_species)

  go_annotations_query <- go_annotations_subject |> 
    dplyr::left_join(ortho_Ath |> dplyr::select(Orthogroup, Gene), by = c("TAIR" = "Gene")) |>
	  dplyr::inner_join(ortho_query |> dplyr::select(Orthogroup, Gene), by = "Orthogroup", relationship = "many-to-many") |>
	  dplyr::distinct(Gene, GO, .keep_all = TRUE) |>
	  dplyr::select(Gene, GO, ONTOLOGY, TERM)

  return(go_annotations_query)
}

run_GO <- function(gene_selected, gene_background, pval_cutoff = 0.05, qval_cutoff = 0.05) {
  gene_selected    <- as.character(gene_selected)
  gene_background  <- as.character(gene_background)

  ego_df <- clusterProfiler::enricher(
    gene = gene_selected,
    universe = gene_background,
    TERM2GENE = go_gene_ontology_table[, c("GO", "Gene")],
    TERM2NAME = go_gene_ontology_table[, c("GO", "TERM")],
    pAdjustMethod = "BH",
    pvalueCutoff = pval_cutoff,
    qvalueCutoff = qval_cutoff,
    minGSSize = 1,
  ) |> tibble::as_tibble() |>
    dplyr::inner_join(go_gene_ontology_table |> dplyr::select(GO, ONTOLOGY) |> distinct(), by = c("ID" = "GO")) |>
    dplyr::arrange(ONTOLOGY, p.adjust)

  ego_df

}


