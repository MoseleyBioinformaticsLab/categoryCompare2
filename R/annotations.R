assign_features <- function(enrichment_obj){
  sig_features <- lapply(enrichment_obj, function(x){
    x@features
  })
  
  all_features <- unique(unlist(sig_features))
  
  sig_matrix <- matrix(FALSE, nrow = length(all_features), ncol = length(sig_features))
  rownames(sig_matrix) <- all_features
  colnames(sig_matrix) <- names(sig_features)
  
  for (isig in names(sig_features)) {
    sig_matrix[sig_features[[isig]], isig] <- TRUE
  }
  
  feature_description <- .annotation_combinations(sig_matrix)@description
  feature_description
}

#' annotation to genes
#' 
#' Creates a tabular output of annotations to genes providing lookup of which
#' genes are contributing to a particular annotation.
#' 
#' @param combined_enrichment combined enrichment object
#' @param annotations which annotations to grab features from
#' @param use_db the annotation database
#' @param input_type what type of gene id was it?
#' @param gene_info what type of info to return for each gene
#' 
#' @export
#' @return data.frame
annotation_gene_table <- function(combined_enrichment,
                                  annotations = NULL,
                                  use_db = NULL,
                                  input_type = "ENTREZID",
                                  gene_info = c("SYMBOL", "GENENAME")){

  annotation_obj <- combined_enrichment@annotation
  if (is.null(annotations)) {
    features <- annotation_obj@annotation_features
  } else {
    annotations <- intersect(annotations, names(annotation_obj@annotation_features))
    features <- annotation_obj@annotation_features[annotations]
  }
  
  
  all_features <- unique(unlist(features))
  feature_info <- suppressMessages(AnnotationDbi::select(use_db, keys = all_features,
                                        keytype = input_type, columns = gene_info))
  
  feature_description <- assign_features(combined_enrichment@enriched)
  
  feature_info <- feature_info[(feature_info[, 1] %in% names(feature_description)), ]
  feature_info$significant <- feature_description[feature_info[, 1]]
  
  out_tables <- lapply(names(features), function(in_annotation){
    header_info <- paste0(in_annotation, " - ", annotation_obj@description[in_annotation])
    annotation_features <- features[[in_annotation]]
    feature_table <- feature_info[(feature_info[, 1] %in% annotation_features), ]
    feature_table <- feature_table[(order(feature_table$significant)), ]
    rownames(feature_table) <- NULL
    list(header = header_info, table = feature_table)
  })
  out_tables
}

#' print table csv
#' 
#' print the annotation gene table to a CSV file
#' 
#' @param annotation_gene_table list of tables
#' @param out_file the file to write to
#' 
#' @export
#' @return NULL
csv_annotation_table <- function(annotation_gene_table, out_file = NULL){
  tmp_file <- tempfile(pattern = "cc2_")
  
  cat(annotation_gene_table[[1]]$header, "\n", file = tmp_file)
  suppressWarnings(write.table(annotation_gene_table[[1]]$table, file = tmp_file, sep = ",",
                   row.names = FALSE, append = TRUE))
  
  for (itable in seq(2, length(annotation_gene_table))) {
    cat("\n\n", annotation_gene_table[[itable]]$header, "\n", file = tmp_file, append = TRUE)
    suppressWarnings(write.table(annotation_gene_table[[itable]]$table, file = tmp_file, sep = ",",
                                 row.names = FALSE, append = TRUE))
    
  }
  if (is.null(out_file)) {
    return(scan(tmp_file, sep = "\n", what = character()))
  } else {
    file.copy(tmp_file, out_file)
    message(paste0("Wrote to ", out_file))
  }
}

#' print table kable
#' 
#' print the annotation gene table in \code{knitr::kable} format
#' 
#' @param annotation_gene_table list of tables
#' @param header_level what header level should the labels be done at?
#' @param cat whether to write it directly, or just return the table for later
#' 
#' @export
#' @return character
kable_annotation_table <- function(annotation_gene_table, header_level = 3, cat = TRUE){
  use_header <- paste0(paste(rep("#", header_level), collapse = ""), " ")
  
  out_table <- character(0)
  
  for (itable in seq(1, length(annotation_gene_table))) {
    out_table <- c(out_table, "\n", "\n")
    out_table <- c(out_table, paste0(use_header, annotation_gene_table[[1]]$header), "\n", "\n")
    out_table <- c(out_table, knitr::kable(annotation_gene_table[[1]]$table), "\n", "\n")
  }
  
  if (cat) {
    cat(out_table, sep = "\n")
  } else {
    return(out_table)
  }
}
