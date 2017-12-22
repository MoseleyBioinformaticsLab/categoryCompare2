#!/usr/bin/Rscript
"
Usage: 
  categoryCompare2.R [--genes=<gene-file>] [--output-directory=<save-location>] [--annotations=<annotation-source>] [--annotation-type=<type>] [--enrichment-test=<enrichment-test>] [--enrichment-direction=<direction>] [--p-adjustment=<p-value adjustment>] [--p-cutoff=<p-value cutoff>] [--count-cutoff=<min-genes>] [--graph=<use-graph>] [--graph-min-edge-weight=<weight-cutoff>] [--graph-communities=<use-communities>] [--graph-vis-type=<visual-node-type>] [--graph-vis-engine=<visualization-engine>]
  categoryCompare2.R [--config=<config-file>]
  categoryCompare2.R [--default-config]
  categoryCompare2.R (-h | --help)
  categoryCompare2.R (-v | --version)

Description: Runs categoryCompare2 on one or more gene lists. Note that there are
a lot of parameters, with several defaults. See the accompanying vignette for a
description of what each one does. The default values can be changed at the command
line, or using a yaml config file.
Options:
  --config=<config-file>                A YAML configuration file [default: NULL]
  --default-config                      Display a default configuration file
  --genes=<gene-file>                   The JSON file containing the genes [default: gene_file.json]
  --geneid-type=<type>                    What type of gene IDs are provided? [default: ENTREZID]
  --output-directory=<save-location>    Where to save the results [default: cc2_results]
  --annotations=<annotation-source>     The annotations to use, either a file or Bioconductor DB [default: org.Hs.eg.db]
  --annotation-type=<type>              If annotations is a Bioconductor DB, then specify the type [default: BP]
  --enrichment-test=<enrichment-test>   What type of test to do [default: hypergeometric]
  --enrichment-direction=<direction>    Do you want over- or under-enrichment [default: over]
  --p-adjustment=<p-value adjustment>   What kind of p-value correction to perform [default: BH]
  --p-cutoff=<p-value cutoff>           What cutoff is required to denote significance? [default: 0.01]
  --count-cutoff=<min-genes>            How many genes need to be annotated to keep the annotation? [default: 2]
  --graph=<use-graph>                   Should the annotation similarity graph be used? [default: TRUE]
  --graph-min-edge-weight=<weight-cutoff> The minimum edge weight to keep graph connections [default: 0.8]
  --graph-communities=<use-communities> Should annotation communities be detected? [default: TRUE]
  --graph-vis-type=<visual-node-type>   How to color the annotation nodes in the graph? [default: pie]
  --graph-vis-engine=<visualization-engine> What software to use to visualize the annotation graph? [default: cytoscape]


" -> doc

default_config <- "genes: gene_file.json
output-directory: cc2_results
annotations: org.Hs.eg.db
annotation-type: BP
enrichment-test: hypergeometric
enrichment-direction: over
p-adjustment: BH
p-cutoff: 0.01
count-cutoff: 2
graph:
  min-edge-weight: 0.8
  communities: TRUE
  vis-type: pie
  vis-engine: cytoscape
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

if (script_options$`default-config`) {
  cat(default_config, sep = "\n")
}

if (script_options$v | script_options$version) {
  message(paste0("This is categoryCompare2 ", packageVersion("categoryCompare2")))
}

#print(script_options)

#script_options <- script_options[!is.null(script_options)]

print(script_options)

# currently, this assumes that *all* the necessary options are set in the config.yml,
# but eventually the goal would be it would just override defaults
if (!(script_options$config %in% "NULL")) {
  tmp_options <- yaml::read_yaml(file = script_options$config)
  print(tmp_options)
  script_options <- tmp_options
  # add code for overriding the defaults
}

# testing: script_options <- yaml::read_yaml(file = "cc2_config.yml")

if (file.exists(script_options$genes)) {
  gene_list <- jsonlite::fromJSON(script_options$genes)
  if (is.null(gene_list$universe)) {
    stop("The gene list does not have a universe / background!")
  } else {
    gene_universe <- gene_list$universe
  }
} else {
  stop("The gene list file ", script_options$genes, " does not exist. Make sure it exists on the search path!")
}

go_types <- c("GO", "BP", "MF", "CC")
go_sub <- c("BP", "MF", "CC")

# check if we are working with an organism db annotation
if (grepl("eg.db", script_options$annotations)) {
  if (!require(script_options$annotations, character.only = TRUE, quietly = TRUE)) {
    stop("The package ", script_options$annotations, " is not installed/available. Try installing it with biocLite('", script_options$annotations, "')")
  } else {
    annotation_src <- eval(parse(text = script_options$annotations))
    
    if (script_options$`annotation-type` %in% go_types) {
      suppressMessages(require("GO.db", character.only = TRUE))
      gene_ann_map <- suppressMessages(AnnotationDbi::select(annotation_src, keys = gene_list$universe,
                                            keytype = script_options$`gene-type`,
                                            columns = "GOALL"))
      
      if (script_options$`annotation-type` %in% go_sub) {
        gene_ann_map <- gene_ann_map[gene_ann_map$ONTOLOGYALL %in% script_options$`annotation-type`, ]
      }
      
      
      ann_gene_list <- split(gene_ann_map[[script_options$`gene-type`]], gene_ann_map[["GOALL"]])
      ann_gene_list <- purrr::map(ann_gene_list, unique)
      ann_description <- suppressMessages(AnnotationDbi::select(GO.db, keys = names(ann_gene_list), columns = "TERM", keytype = "GOID")$TERM)
      names(ann_description) <- names(ann_gene_list)
      
      if (script_options$`annotation-type` %in% "GO") {
        go_ontology_map <- unique(gene_ann_map[, c("GOALL", "ONTOLOGYALL")])
        go_ontology <- go_ontology_map$ONTOLOGYALL
        names(go_ontology) <- go_ontology_map$GOALL
        go_ontology <- go_ontology[names(ann_description)]
        ann_description <- paste0(go_ontology, ":", ann_description)
      }
      
      
      annotation_obj <- categoryCompare2::annotation(annotation_features = ann_gene_list,
                                                     description = ann_description,
                                                     type = script_options$`annotation-type`)
    }
  }
}

# now that we are done making annotations, can remove the universe
gene_list$universe <- NULL

gene_enrichments <- purrr::map(gene_list, function(in_genes){
  hypergeometric_feature_enrichment(
    new("hypergeom_features", significant = in_genes,
        universe = gene_universe, annotation = annotation_obj),
    p_adjust = script_options$`p-adjustment`
  )
})

combined_enrichments <- combine_enrichments(gene_enrichments)

if (!(script_options$`p-adjustment` %in% "none")) {
  p_cutoff_column <- "padjust"
} else {
  p_cutoff_column <- "p"
}

p_cutoff_value <- script_options$`p-cutoff`
p_cutoff_direction <- "<="

count_cutoff_column <- "counts"
count_cutoff_value <- script_options$`count-cutoff`
count_cutoff_direction <- ">="

count_call_info <- list(fun = count_cutoff_direction, var_1 = count_cutoff_column, var_2 = count_cutoff_value)
p_call_info <- list(fun = p_cutoff_direction, var_1 = p_cutoff_column, var_2 = p_cutoff_value)

significant_calls <- list(counts = count_call_info, pvalues = p_call_info)

combined_enrichments <- combined_significant_calls(combined_enrichments, significant_calls)

