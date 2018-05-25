#!/usr/bin/Rscript
"
Usage: 
  filter_and_group.R --enrichment-directory=<enrichment_directory> [--p-cutoff=<max-p-value>] [--adjusted-p-values=<use-adjusted-p-values>] [--count-cutoff=<min-genes>] [--group=<do-grouping>] [--similarity-cutoff=<minimum similarity>] [--grouping-algorithm=<group-algorithm>] [--table-file=<table-file>]
  filter_and_group.R (-h | --help)

Description: Given an enrichment directory, will load up the enrichment results, and then performs filtering
of results to denote significant annotations in each group, and optionally tries to group annotations based
on similarity and a grouping algorithm. Results are written to the table file.

Valid choices for the grouping algorithm include:
  edge_betweenness
  fast_greedy
  label_prop
  leading_eigen
  louvain
  optimal
  spinglass
  walktrap

Options:
  --enrichment-directory=<enrichment_directory>   where the enrichment results are found
  --p-cutoff=<max-p-value>                        the maximum p-value to consider signifcant [default: 0.01] 
  --adjusted-p-values=<use-adjusted-p-values>     should adjusted p-values be used if they exist? [default: TRUE]
  --count-cutoff=<min-genes>                      minimum number of genes annotated [default: 2] 
  --group=<do-grouping>                           should grouping of annotations be attempted [default: TRUE]
  --similarity-cutoff=<minimum similarity>        minimum similarity measure to consider annotations linked [default: 0] 
  --grouping-algorithm=<group-algorithm>          what algorithm should be used to find the groups [default: walktrap]
  --table-file=<table-file>                       the results file to save the results
  --network-file=<network-file>                   if desired, save the network as well [default: NULL ]

" -> doc

library(methods)
library(docopt)
library(tools)
library(jsonlite)
library(yaml)
suppressMessages(library(categoryCompare2))

script_options <- docopt(doc)

main <- function(script_options){
  print(script_options)
  
  grouping_algorithms <- c("walktrap" = "cluster_walktrap",
                           "spinglass" = "cluster_spinglass",
                           "optimal" = "cluster_optimal",
                           "louvain" = "cluster_louvain",
                           "leading_eigen" = "cluster_leading_eigen",
                           "label_prop" = "cluster_label_prop",
                           "fast_greedy" = "cluster_fast_greedy",
                           "edge_betweenness" = "cluster_edge_betweenness")
  
  # fail before doing enrichment!
  p_cutoff_value <- as.double(script_options$`p-cutoff`)
  if (is.na(p_cutoff_value)) {
    stop("The p-cutoff MUST be a number, or Inf")
  }
  p_cutoff_direction <- "<="
  
  count_cutoff_column <- "counts"
  count_cutoff_value <- as.double(script_options$`count-cutoff`)
  if (is.na(count_cutoff_value)) {
    stop("The count-cutoff MUST be a number, or Inf")
  }
  count_cutoff_direction <- ">="
  
  if (script_options$group) {
    similarity_cutoff <- as.double(script_options$`similarity-cutoff`)
    if (is.na(similarity_cutoff)) {
      stop("The similarity-cutoff MUST be a number!")
    }
  }
  
  
  if (!(script_options$`adjusted-p-values` %in% "FALSE")) {
    p_cutoff_column <- "padjust"
  } else {
    p_cutoff_column <- "p"
  }
  
  
  count_call_info <- list(fun = count_cutoff_direction, var_1 = count_cutoff_column, var_2 = count_cutoff_value)
  p_call_info <- list(fun = p_cutoff_direction, var_1 = p_cutoff_column, var_2 = p_cutoff_value)
  
  significant_calls <- list(counts = count_call_info, pvalues = p_call_info)
  
  enrichment_file <- file.path(script_options$`enrichment-directory`, "enrichments.rds")
  if (file.exists(enrichment_file)) {
    enrichments <- readRDS(enrichment_file)
  } else {
    stop("Enrichment results file does not exist!")
  }
  
  significant <- combined_significant_calls(enrichments, significant_calls)
  message("Significant Annotations:")
  print(significant@statistics@significant)
  
  table_file <- file.path(script_options$`enrichment-directory`, script_options$`table-file`)
  
  if (!as.logical(script_options$group)) {
    results_table <- generate_table(significant)
    write.table(results_table, file = table_file, sep = "\t", row.names = FALSE, col.names = TRUE)
  } else {
    similarity_graph <- generate_annotation_graph(significant)
    
    if (similarity_cutoff > 0) {
      similarity_graph <- remove_edges(similarity_graph, similarity_cutoff)
    }
    similarity_graph
    
    significant_assignments <- annotation_combinations(similarity_graph)
    graph_communities <- assign_communities(similarity_graph)
    community_labels <- label_communities(graph_communities, enrichments@annotation)
    
    results_table <- table_from_graph(similarity_graph, significant_assignments, community_labels)
    write.table(results_table, file = table_file, sep = "\t", row.names = FALSE, col.names = TRUE)
  }
}

main(script_options)
