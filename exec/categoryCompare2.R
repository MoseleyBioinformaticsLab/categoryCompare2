#!/usr/bin/Rscript
"
Usage: 
  categoryCompare2.R [--features=<feature-file>] [--output-directory=<save-location>] [--annotations=<annotation-source>] [--annotation-type=<type>] [--enrichment-test=<enrichment-test>] [--enrichment-direction=<direction>] [--p-adjustment=<p-value adjustment>] [--p-cutoff=<p-value cutoff>] [--count-cutoff=<min-genes>] [--graph=<use-graph>] [--graph-min-edge-weight=<weight-cutoff>] [--graph-communities=<use-communities>] [--graph-vis-type=<visual-node-type>] [--graph-vis-engine=<visualization-engine>]
  categoryCompare2.R [--config=<config-file>]
  categoryCompare2.R [--default-config]
  categoryCompare2.R (-h | --help)
  categoryCompare2.R (-v | --version)

Description: Runs categoryCompare2 on one or more feature lists. Note that there are
a lot of parameters, with several defaults. See the accompanying vignette for a
description of what each one does. The default values can be changed at the command
line, or using a yaml config file.
Options:
  --config=<config-file>                A YAML configuration file [default: NULL]
  --default-config                      Display a default configuration file
  --features=<feature-file>             The JSON file containing the features (genes) [default: features.json]
  --output-directory=<save-location>    Where to save the results [default: cc2_results]
  --annotations=<annotation-source>     The annotations to use, as a file [default: annotations.json]
  --enrichment-test=<enrichment-test>   What type of test to do [default: hypergeometric]
  --enrichment-direction=<direction>    Do you want over- or under-enrichment [default: over]
  --p-adjustment=<p-value adjustment>   What kind of p-value correction to perform [default: BH]
  --p-cutoff=<p-value cutoff>           What cutoff is required to denote significance (numeric)? [default: 0.001]
  --count-cutoff=<min-genes>            How many genes need to be annotated to keep the annotation (numeric)?  [default: 2]

" -> doc

default_config <- "features: feature_file.json
output-directory: cc2_results
annotations: annotations.json
enrichment-test: hypergeometric
enrichment-direction: over
p-adjustment: BH
p-cutoff: 0.001
count-cutoff: 2
"

# testing:
# setwd("executable_related")
# script_options <- list(file1 = "test_data/10_entrez.txt", file2 = "test_data/48_entrez.txt", universe = "test_data/universe_entrez.txt")
# 

library(methods)
library(docopt)
library(tools)
library(jsonlite)
library(yaml)
suppressMessages(library(categoryCompare2))


script_options <- docopt(doc)

main <- function(script_options){
  
  #print(script_options)
  
  #browser(expr = TRUE)
  if (!is.null(script_options$`default-config`)) {
    if (script_options$`default-config`) {
      cat(default_config, sep = "\n")
      return()
    }
    
  }
  
  
  
  if (!is.null(script_options$version)) {
    if (script_options$version) {
      message(paste0("This is categoryCompare2 ", packageVersion("categoryCompare2")))
      return()
    }
    
  }
  
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
  
  if (!(script_options$`p-adjustment` %in% "none")) {
    p_cutoff_column <- "padjust"
  } else {
    p_cutoff_column <- "p"
  }
  
  
  count_call_info <- list(fun = count_cutoff_direction, var_1 = count_cutoff_column, var_2 = count_cutoff_value)
  p_call_info <- list(fun = p_cutoff_direction, var_1 = p_cutoff_column, var_2 = p_cutoff_value)
  
  significant_calls <- list(counts = count_call_info, pvalues = p_call_info)
  
  #print(script_options)
  
  #script_options <- script_options[!is.null(script_options)]
  
  
  
  # currently, this assumes that *all* the necessary options are set in the config.yml,
  # but eventually the goal would be it would just override defaults
  if (!is.null(script_options$config)) {
    if (!(script_options$config %in% "NULL")) {
      tmp_options <- yaml::read_yaml(file = script_options$config)
      print(tmp_options)
      script_options <- tmp_options
    }
    # add code for overriding the defaults
  }
  
  # testing: script_options <- yaml::read_yaml(file = "cc2_config.yml")
  
  if (file.exists(script_options$features)) {
    feature_list <- jsonlite::fromJSON(script_options$features)
    if (is.null(feature_list$universe)) {
      stop("The feature list does not have a universe / background!")
    } else {
      feature_universe <- feature_list$universe
      feature_list$universe <- NULL
    }
  } else {
    stop("The feature list file ", script_options$features, " does not exist. Make sure it exists on the search path!")
  }
  
  # read in the annotations and create the annotation object
  annotation_obj <- json_2_annotation(script_options$annotations)
  
  annotation_features <- unique(unlist(annotation_obj@annotation_features))
  features_in_annotations <- vapply(feature_list, function(x){
    sum(x %in% annotation_features)
  }, numeric(1))
  
  if (any(vapply(features_in_annotations, function(x){sum(x) == 0}, logical(1)))) {
    stop("One or more of your feature lists has 0 matches to the features in the annotations!")
  }
  
  gene_enrichments <- lapply(feature_list, function(in_genes){
    hypergeometric_feature_enrichment(
      new("hypergeom_features", significant = in_genes,
          universe = feature_universe, annotation = annotation_obj),
      p_adjust = script_options$`p-adjustment`
    )
  })
  
  combined_enrichments <- combine_enrichments(gene_enrichments)
  
  
  combined_significant <- combined_significant_calls(combined_enrichments, significant_calls)
  
  
  
  results_table <- generate_table(combined_significant)
  
  if (!dir.exists(script_options$`output-directory`)) {
    dir.create(script_options$`output-directory`, recursive = TRUE)
  }
  #saveRDS(combined_significant, file = file.path(script_options$`output-directory`, "combined_significant.rds"))
  #print(combined_significant)
  write.table(results_table, file = file.path(script_options$`output-directory`, "full_table.txt"),
              sep = "\t", row.names = FALSE, col.names = TRUE)
  
  # next up will be code useful for handling wanting to make a graph
  # significant_graph <- generate_annotation_graph(combined_significant)
  # significant_graph <- remove_edges(significant_graph, script_options$graph$`min-edge-weight`)
  
}

# testing: script_options <- yaml::read_yaml(file = "cc2_config.yml")
main(script_options)
