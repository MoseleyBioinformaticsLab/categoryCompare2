#!/usr/bin/env Rscript
"
Usage: 
  run_enrichment.R [--features=<feature-file>] [--annotations=<annotation-source>] [--enrichment-test=<enrichment-test>] [--enrichment-direction=<direction>] [--p-adjustment=<p-value adjustment>] [--output-file=<output_file>] [--text-only=<text-only>]
  run_enrichment.R [--config=<config-file>]
  run_enrichment.R [--default-config]
  run_enrichment.R (--version | -v)
  run_enrichment.R (--help | -h)

Description: Runs categoryCompare2 on one or more feature lists. Note that there are
a lot of parameters, with several defaults. See the accompanying vignette for a
description of what each one does. The default values can be changed at the command
line, or using a yaml config file.

Unless --output-text-only=TRUE, then both a tab delimited text file and a binary
`rds` file will be written. The `rds` file is necessary for `filter_and_group.R`.
If you don't want to do filtering and grouping, then you probably want to use
--output-text-only=TRUE.

Options:
  --config=<config-file>                A YAML configuration file [default: NULL]
  --default-config                      Display a default configuration file
  --features=<feature-file>             The JSON file containing the features (genes) [default: features.json]
  --annotations=<annotation-source>     The annotations to use, as a file [default: annotations.json]
  --enrichment-test=<enrichment-test>   What type of test to do [default: hypergeometric]
  --enrichment-direction=<direction>    Do you want over- or under-enrichment [default: over]
  --p-adjustment=<p-value adjustment>   What kind of p-value correction to perform [default: BH]
  --output-file=<save-location>         Where to save the results [default: cc2_results.txt]
  --text-only=<text-only>               Should only the text file be generated? [default: FALSE]


" -> doc

default_config <- "features: features.json
output-file: cc2_results/cc2_results.txt
annotations: annotations.json
enrichment-test: hypergeometric
enrichment-direction: over
p-adjustment: BH
"

# testing:
# setwd("executable_related")
# script_options <- list(file1 = "test_data/10_entrez.txt", file2 = "test_data/48_entrez.txt", universe = "test_data/universe_entrez.txt")
#
# setwd("executable_related")
# script_options <- list(features = "test_data/gsea_features.json", output_file = "test_data/cc2_gsea_results.txt", enrichment_test = "gsea",
# annotations = "test_data/go_entrez_annotations.json", output_file = "test_data/gsea_output.txt")

library(methods)
library(docopt)
library(tools)
library(jsonlite)
library(yaml)
suppressMessages(library(categoryCompare2))


script_options <- docopt(doc)

main <- function(script_options) {
  #print(script_options)

  #browser(expr = TRUE)
  if (!is.null(script_options$default_config)) {
    if (script_options$default_config) {
      cat(default_config, sep = "\n")
      return()
    }
  }

  if (!is.null(script_options$version)) {
    if (script_options$version) {
      message(paste0(
        "This is categoryCompare2 ",
        packageVersion("categoryCompare2")
      ))
      return()
    }
  }

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
  } else {
    stop(
      "The feature list file ",
      script_options$features,
      " does not exist. Make sure it exists on the search path!"
    )
  }

  #print(script_options$output_file)
  output_dir <- dirname(script_options$output_file)

  if (!dir.exists(output_dir)) {
    message(paste0("Creating directory ", output_dir))
    dir.create(output_dir, recursive = TRUE)
  }

  # read in the annotations and create the annotation object
  annotation_obj <- json_2_annotation(script_options$annotations)

  run_hypergeometric <- function(script_options, feature_list, annotation_obj) {
    if (is.null(feature_list$universe)) {
      stop(
        "The feature list does not have a universe / background, which is required for hypergeometric enrichment."
      )
    } else {
      feature_universe <- feature_list$universe
      feature_list$universe <- NULL
    }
    annotation_features <- unique(unlist(annotation_obj@annotation_features))
    features_in_annotations <- vapply(
      feature_list,
      function(x) {
        sum(x %in% annotation_features)
      },
      numeric(1)
    )

    if (
      any(vapply(
        features_in_annotations,
        function(x) {
          sum(x) == 0
        },
        logical(1)
      ))
    ) {
      stop(
        "One or more of your feature lists has 0 matches to the features in the annotations!"
      )
    }

    gene_enrichments <- lapply(feature_list, function(in_genes) {
      hypergeometric_feature_enrichment(
        new(
          "hypergeometric_features",
          significant = in_genes,
          universe = feature_universe,
          annotation = annotation_obj
        ),
        p_adjust = script_options$p_adjustment
      )
    })
    return(gene_enrichments)
  }

  run_gsea <- function(script_options, feature_list, annotation_obj) {
    gene_enrichments <- lapply(feature_list, function(in_genes) {
      rank_values <- in_genes$rank
      names(rank_values) <- in_genes$feature
      gsea_feature_enrichment(
        new("gsea_features", ranks = rank_values, annotation = annotation_obj)
      )
    })
  }

  if (script_options$enrichment_test %in% "hypergeometric") {
    gene_enrichments <- run_hypergeometric(
      script_options,
      feature_list,
      annotation_obj
    )
  } else if (script_options$enrichment_test %in% "gsea") {
    gene_enrichments <- run_gsea(script_options, feature_list, annotation_obj)
  }

  combined_enrichments <- combine_enrichments(gene_enrichments)

  rds_file <- paste0(
    tools::file_path_sans_ext(script_options$output_file),
    ".rds"
  )

  if (!as.logical(script_options$text_only)) {
    saveRDS(combined_enrichments, file = file.path(rds_file))
  }

  results_table <- generate_table(combined_enrichments)

  #saveRDS(combined_significant, file = file.path(script_options$`output-directory`, "combined_significant.rds"))
  #print(combined_significant)
  has_leading <- grepl("leading_edge", colnames(results_table))
  if (any(has_leading)) {
    results_table = results_table[, !has_leading]
  }
  write.table(
    results_table,
    file = file.path(script_options$output_file),
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

  # next up will be code useful for handling wanting to make a graph
  # significant_graph <- generate_annotation_graph(combined_significant)
  # significant_graph <- remove_edges(significant_graph, script_options$graph$`min-edge-weight`)
}

# testing: script_options <- yaml::read_yaml(file = "cc2_config.yml")
main(script_options)
