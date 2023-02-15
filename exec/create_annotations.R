#!/usr/bin/Rscript
"
Usage: 
  create_annotations.R [--input=<input-file>] [--orgdb=<organism_db>] [--feature-type=<feature-type>] [--annotation-type=<annotation-type>] [--json=<jsonfile>] 
  create_annotations.R (-h | --help)
  create_annotations.R (-v | --version)

Description: Either take a feature to annotation JSON file, and reverse them
  to create the needed annotation object for categoryCompare, or generate the
  annotations using a Bioconductor organism or chip database.

Options:
  --input=<input-file>                  An input JSON file of features to annotations (optional)
  --orgdb=<organism_db>                 A Bioconductor organism or chip database to use (optional)
  --feature-type=<feature-type>         What type of features to use from the database [default: ENTREZID]
  --annotation-type=<annotation-type>   The annotation type [default: GO]
  --json=<jsonfile>                     The output file [default: annotations.json]
  
" -> doc

library(methods)
library(docopt)
library(tools)
library(jsonlite)
suppressMessages(library(categoryCompare2))
script_options <- docopt(doc)

main <- function(script_options){
  null_input <- is.null(script_options$input)
  null_orgdb <- is.null(script_options$orgdb)
  
  if (!is.null(script_options$version)) {
    if (script_options$version) {
      message(paste0("This is categoryCompare2 ", packageVersion("categoryCompare2")))
      return()
    }
    
  }
  
  if (null_input && null_orgdb) {
    stop("One of --input or --orgdb have to be defined!")
  }
  
  
  
  if (!is.null(script_options$input)) {
    json_annotation_reversal(script_options$input, out_file = script_options$json, 
                        feature_type = script_options$feature_type,
                        annotation_type = script_options$annotation_type)
  }
  
  if (!is.null(script_options$orgdb)) {
    annotation_obj <- get_db_annotation(orgdb = script_options$orgdb, feature_type = script_options$feature_type,
                                        annotation_type = script_options$annotation_type)
    
    annotation_2_json(annotation_obj, script_options$json)
  }
}

main(script_options)
