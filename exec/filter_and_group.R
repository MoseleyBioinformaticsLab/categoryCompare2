#!/usr/bin/env Rscript
"
Usage: 
  filter_and_group.R [--enrichment-results=<enrichment_results>] [--p-cutoff=<max-p-value>] [--adjusted-p-values=<use-adjusted-p-values>] [--count-cutoff=<min-features>] [--similarity-file=<similarity-file>] [--similarity-cutoff=<minimum similarity>] [--grouping-algorithm=<group-algorithm>] [--table-file=<table-file>]
  filter_and_group.R (-h | --help)
  filter_and_group.R (-v | --version)

Description: Given an enrichment directory, will load up the enrichment results, and then performs filtering
of results to denote significant annotations in each group, and optionally tries to group annotations based
on similarity and a grouping algorithm. Results are written to the table file.

Annotation similarities are calculated and annotations grouped if a file is specified.
If the file exists, then the annotation similarities will not be re-calculated, but they will
be loaded from the file and used directly.

Valid choices for the grouping algorithm include:
  edge_betweenness
  fast_greedy
  label_prop
  leading_eigen
  louvain
  optimal
  spinglass
  walktrap

Note that although the enrichment-results is specified as a `.txt`, this script
will look for a matching `.rds` file, and will use that instead. All of the information
needed for filtering and grouping is not captured in the text file representation.

Annotation similarities are calculated and group

Options:
  --enrichment-results=<enrichment_results>       where the enrichment results are found [default: cc2_results.txt]
  --p-cutoff=<max-p-value>                        the maximum p-value to consider significant [default: 0.01] 
  --adjusted-p-values=<use-adjusted-p-values>     should adjusted p-values be used if they exist? [default: TRUE]
  --count-cutoff=<min-features>                   minimum number of features annotated [default: 2] 
  --similarity-file=<similarity-file>             where should annotation similarities be saved? [default: annotation_similarity.rds]
  --similarity-cutoff=<minimum similarity>        minimum similarity measure to consider annotations linked [default: 0] 
  --grouping-algorithm=<group-algorithm>          what algorithm should be used to find the groups [default: walktrap]
  --table-file=<table-file>                       the results file to save the results [default: cc2_results_grouped.txt]
  --network-file=<network-file>                   if desired, save the network as well [default: NULL]

" -> doc

library(methods)
library(docopt)
library(tools)
library(jsonlite)
library(yaml)
suppressMessages(library(categoryCompare2))

script_options <- docopt(doc)

main <- function(script_options) {
  #print(script_options)

  grouping_algorithms <- c(
    "walktrap" = "cluster_walktrap",
    "spinglass" = "cluster_spinglass",
    "optimal" = "cluster_optimal",
    "louvain" = "cluster_louvain",
    "leading_eigen" = "cluster_leading_eigen",
    "label_prop" = "cluster_label_prop",
    "fast_greedy" = "cluster_fast_greedy",
    "edge_betweenness" = "cluster_edge_betweenness"
  )

  if (!is.null(script_options$version)) {
    if (script_options$version) {
      message(paste0(
        "This is categoryCompare2 ",
        packageVersion("categoryCompare2")
      ))
      return()
    }
  }

  # fail before doing filtering!
  enrichment_file <- paste0(
    tools::file_path_sans_ext(script_options$enrichment_results),
    ".rds"
  )
  if (file.exists(enrichment_file)) {
    enrichments <- readRDS(enrichment_file)
  } else {
    stop("Enrichment results RDS file does not exist!")
  }

  method_type <- enrichments@enriched[[1]]@statistics@method

  p_cutoff_value <- as.double(script_options$p_cutoff)
  if (is.na(p_cutoff_value)) {
    stop("The p-cutoff MUST be a number, or Inf")
  }
  p_cutoff_direction <- "<="

  count_cutoff_column <- "counts"
  count_cutoff_value <- as.double(script_options$count_cutoff)
  if (is.na(count_cutoff_value)) {
    stop("The count-cutoff MUST be a number, or Inf")
  }
  count_cutoff_direction <- ">="

  if (!is.null(script_options$similarity_file)) {
    similarity_cutoff <- as.double(script_options$similarity_cutoff)
    if (is.na(similarity_cutoff)) {
      stop("The similarity-cutoff MUST be a number!")
    }
  }

  if (!(script_options$adjusted_p_values %in% "FALSE")) {
    p_cutoff_column <- "padjust"
  } else {
    p_cutoff_column <- "p"
  }

  count_call_info <- list(
    fun = count_cutoff_direction,
    var_1 = count_cutoff_column,
    var_2 = count_cutoff_value
  )
  p_call_info <- list(
    fun = p_cutoff_direction,
    var_1 = p_cutoff_column,
    var_2 = p_cutoff_value
  )

  if (method_type %in% "hypergeometric") {
    significant_calls <- list(counts = count_call_info, pvalues = p_call_info)
  } else if (method_type %in% "gsea") {
    significant_calls <- list(pvalues = p_call_info)
  }

  significant <- combined_significant_calls(enrichments, significant_calls)
  message("Significant Annotations:")
  print(significant@statistics@significant)

  n_significant = unlist(significant@statistics@significant@significant) |>
    sum()

  if (n_significant == 0) {
    message(
      "Nothing significant in any enrichments, stopping.\nConsider adjusting `--p-cutoff` or `--count-cutoff`."
    )
    return()
  }

  table_dir <- dirname(script_options$table_file)
  if (!dir.exists(table_dir)) {
    message(paste0("Creating directory: ", table_dir))
    dir.create(table_dir, recursive = TRUE)
  }

  if (script_options$similarity_file == "NULL") {
    results_table <- generate_table(significant)
    write.table(
      results_table,
      file = script_options$table_file,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
  } else {
    if (file.exists(script_options$similarity_file)) {
      message(paste0(
        "Loading annotation similarities from: ",
        script_options$similarity_file
      ))
      similarity_graph <- readRDS(script_options$similarity_file)
    } else {
      similarity_graph <- generate_annotation_graph(significant)
      similarity_dir <- dirname(script_options$similarity_file)
      if (!dir.exists(similarity_dir)) {
        message(paste0("Creating directory: ", similarity_dir))
        dir.create(similarity_dir, recursive = TRUE)
      }
      message(paste0(
        "Saving annotation similarities in: ",
        script_options$similarity_file
      ))
      saveRDS(similarity_graph, file = script_options$similarity_file)
    }

    if (similarity_cutoff > 0) {
      similarity_graph <- remove_edges(similarity_graph, similarity_cutoff)
    }
    similarity_graph

    significant_assignments <- annotation_combinations(similarity_graph)
    graph_communities <- assign_communities(similarity_graph)
    community_labels <- label_communities(
      graph_communities,
      enrichments@annotation
    )

    results_table <- table_from_graph(
      similarity_graph,
      significant_assignments,
      community_labels
    )
    write.table(
      results_table,
      file = script_options$table_file,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE
    )
  }
}

main(script_options)
