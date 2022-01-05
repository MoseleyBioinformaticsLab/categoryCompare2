#' get significant annotations calls
#' 
#' In the case where we have a \code{\link{combined_enrichment}} and we want
#' to get all of the significant annotations from each of them, and put them
#' together so we can start doing real meta-analysis.
#' 
#' Note that this function returns the original \code{\link{combined_enrichment}} object with a modified
#' \code{\link{combined_statistics}} slot where the significant annotations have been added in. 
#' 
#' @param in_results a \code{\link{combined_enrichment}} object
#' @param queries a list of queries that can form a call object
#' 
#' @return \code{\link{combined_enrichment}} object
#' @export
combined_significant_calls <- function(in_results, queries){
  all_measured <- lapply(in_results@enriched,
                         function(x){x@statistics@annotation_id})
  
  all_significant <- lapply(in_results@enriched,
                            function(x){get_significant_annotations_calls(x@statistics, queries)})
  
  annotation_measured <- unique(unlist(all_measured))
  n_measured <- length(annotation_measured)
  n_enriched <- length(in_results@enriched)
  
  out_measured <- matrix(FALSE, n_measured, n_enriched)
  rownames(out_measured) <- annotation_measured
  colnames(out_measured) <- names(all_measured)
  
  out_significant <- out_measured
  
  for (i_meas in names(all_measured)){
    out_measured[all_measured[[i_meas]], i_meas] <- TRUE
  }
  
  for (i_meas in names(all_significant)){
    out_significant[all_significant[[i_meas]], i_meas] <- TRUE
  }
  
  sig_annotation <- new("significant_annotations",
                        significant = out_significant,
                        measured = out_measured,
                        sig_calls = sapply(queries, function(x){
                          paste0(x$var_1, " ", x$fun, " ", x$var_2)
                        }))
  
  in_results@statistics@significant <- sig_annotation
  
  in_results
}

#' get significant annotations calls
#' 
#' In the case where we have a \code{\link{statistical_results}} and we want
#' to get all of the significant annotations from it
#' 
#' @param in_results a \code{\link{statistical_results}} object
#' @param queries a list of queries that can form a call object
#' 
#' @return vector of significant annotation_id's
#' @export
get_significant_annotations_calls <- function(in_results, queries){
  out_ids <- in_results@annotation_id
  
  sig_entries <- multi_query_list_calls(in_results@statistic_data, queries)
  
  out_ids[sig_entries]
  
}

multi_query_list_calls <- function(list_data, queries){
  n_query <- length(queries)
  
  query_vars <- purrr::map_chr(queries, "var_1")
  if (sum(query_vars %in% names(list_data)) != length(query_vars)) {
    stop("Query objects are not all present in your data!")
  }
  
  query_results <- lapply(queries, function(in_query){
    use_call <- call(in_query$fun, list_data[[in_query$var_1]], in_query$var_2)
    eval(use_call)
  })
  
  # how many objects do we have in each query
  # they should all be the same to allow merging
  n_objects <- unique(unlist(lapply(query_results, length)))
  
  if (length(n_objects) != 1){
    stop("cannot merge queries on objects of different length", call. = FALSE)
  }
  
  result_logical <- rep(TRUE, n_objects)
  
  for (i_query in seq(1, n_query)){
    result_logical <- result_logical & query_results[[i_query]]
  }
  
  names(result_logical) <- NULL
  result_logical
}
