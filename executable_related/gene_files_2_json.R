#!/usr/bin/Rscript
"
Usage: 
  gene_files_2_json.R [--json=<jsonfile>] [--file1=<file1>] [--file2=<file2>] [--file3=<file3>] [--file4=<file4>] [--universe=<universe>]
  gene_files_2_json.R (-h | --help)

Description: Concatenates multiple gene list files into a single JSON file. Each
gene list will be named according to the file name that it came from. So, for example
if you did: 

./gene_files_2_json.R --file1=treatment1.txt --file2=treatment2.txt --universe=universe.txt

Then there would be *treatment1*, *treatment2*, and *universe* in the JSON file

Note that if *universe* is not supplied, then it will be the combination of all
of the other gene lists supplied.

Options:
  --json=<jsonfile>     The JSON file to save to [default: gene_file.json]
  --file1=<file1>       The first list of genes (optional)
  --file2=<file2>       The second list of genes (optional)
  --file3=<file3>       The third list of genes (optional)
  --file4=<file4>       The fourth list of genes (optional)
  --universe=<universe> All the genes measured (optional)
  
" -> doc

# testing:
# setwd("executable_related")
# script_options <- list(file1 = "test_data/10_entrez.txt", file2 = "test_data/48_entrez.txt", universe = "test_data/universe_entrez.txt")
# 

library(methods)
library(docopt)
library(tools)
library(jsonlite)
script_options <- docopt(doc)

#print(script_options)

script_options <- script_options[!is.null(script_options)]
file_args_loc <- grepl("^file", names(script_options)) | grepl("^universe", names(script_options))
file_args <- script_options[file_args_loc]

get_gene_lists <- function(file_list){
  file_not_universe <- unlist(file_list[!(names(file_list) %in% "universe")])
  
  condition_names <- basename(file_not_universe)
  condition_names <- gsub(paste0(".", file_ext(condition_names[1])), "", condition_names)
  
  file_data <- lapply(file_not_universe, function(x){
    readLines(x)
  })
  names(file_data) <- condition_names
  
  if (is.null(file_list$universe)) {
    file_data$universe <- unique(unlist(file_data))
  } else {
    file_data$universe <- readLines(file_list$universe)
  }
  
  file_data
}

gene_lists <- get_gene_lists(file_args)

cat(jsonlite::toJSON(gene_lists, pretty = TRUE), file = script_options$json)
