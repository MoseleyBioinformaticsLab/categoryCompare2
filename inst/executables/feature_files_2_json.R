#!/usr/bin/Rscript
"
Usage: 
  feature_files_2_json.R [--json=<jsonfile>] [--file1=<file1>] [--file2=<file2>] [--file3=<file3>] [--file4=<file4>] [--universe=<universe>]
  feature_files_2_json.R (-h | --help)

Description: Concatenates multiple feature list files into a single JSON file. Each
feature list will be named according to the file name that it came from. So, for example
if you did: 

feature_files_2_json.R --file1=treatment1.txt --file2=treatment2.txt --universe=universe.txt

Then there would be *treatment1*, *treatment2*, and *universe* in the JSON file

Note that if *universe* is not supplied, then it will be the combination of all
of the other feature lists supplied.

Options:
  --json=<jsonfile>     The JSON file to save to [default: features.json]
  --file1=<file1>       The first list of features (optional)
  --file2=<file2>       The second list of features (optional)
  --file3=<file3>       The third list of features (optional)
  --file4=<file4>       The fourth list of features (optional)
  --universe=<universe> All the features measured (optional)
  
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


check_feature_lists <- function(feature_list, name_list, universe){
  n_not <- length(setdiff(feature_list, universe))
  
  if (n_not > 0) {
    warn_message <- paste0("There are ", n_not, " features from ", name_list, " list not in the Universe list!")
    warning(warn_message)
  }
}


get_feature_lists <- function(file_list){
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
  
  not_universe <- file_data[!(names(file_data) %in% "universe")]
  invisible(purrr::imap(not_universe, function(.x, .y){
    check_feature_lists(.x, .y, file_data$universe)
  }))
  
  file_data
}

feature_lists <- get_feature_lists(file_args)

cat(jsonlite::toJSON(feature_lists, pretty = TRUE), file = script_options$json)
