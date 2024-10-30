#' GSEA feature class
#' 
#' class to hold features undergoing GSEA
#' 
#' @slot ranks a named vector of ranks
#' @slot annotation annotation object
#' 
#' @export
setClass("gsea_features",
         slots = list(ranks = "ANY",
                      annotation = "annotation"))

#' do GSEA
#' 
#' Performs gene-set enrichment analysis using the `fgsea` package.
#' 
#' @param gsea_features a GSEA features object
#' @param min_features the minimum number of features for an annotation (default = 15)
#' @param max_features the maximum number of features for an annotation (default = 500)
#' @param return_type  what type of object should be returned? ("cc2" or "fgsea")
#' 
#' @details The runtime is dependent on the maximum size of the provided annotation,
#'   so the authors of `fgsea` recommend a maximum size of 500. In addition, to calculate
#'   statistics, a minimum size of annotated features are required. Going below 15 may
#'   not be advised. If you want to use other `fgsea` functions, it is recommended to
#'   set `return_type = "fgsea"`. Otherwise, you should keep the default of "cc2".
#' 
#' @seealso fgsea::fgsea
#' 
#' @export
#' @return enriched_result
gsea_feature_enrichment = function(gsea_features, 
                                   min_features = 15, 
                                   max_features = 500,
                                   return_type = "cc2",
                                   ...)
{
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("The 'fgsea' package must be installed!")
  }

  gsea_res = fgsea::fgsea(pathways = gsea_features@annotation@annotation_features,
                          stats = gsea_features@ranks,
                          minSize = min_features,
                          maxSize = max_features,
                          ...)
  
  if (return_type %in% "fgsea") {
    return(gsea_res)
  } else if (return_type %in% "cc2") {
    annotation_id = gsea_res$pathway
    gsea_res$pathway = NULL
    rename_list = c(p = "pval", padjust = "padj", log2err = "log2err", es = "ES", nes = "NES", size = "size", leading_edge = "leadingEdge")
    gsea_res = as.list(dplyr::rename(gsea_res, dplyr::all_of(rename_list)))
    gsea_res = lapply(gsea_res, function(x){
      names(x) = annotation_id
      x
    })
    stats_results = new("statistical_results",
                        statistic_data = gsea_res,
                        annotation_id = annotation_id,
                        method = "gsea")
    out_enrich = new("enriched_result",
                      ranks = gsea_features@ranks,
                      statistics = stats_results,
                      annotation = gsea_features@annotation)
    return(out_enrich)
  } else {
    stop("return_type not recognized.")
  }
    
}

#' convert enriched object
#' 
#' Takes an `enriched_result`, and converts it to the table expected by `fgsea`. 
#' This should only be done on those that have `gsea` as the *Enrichment Method*.
#' 
#' @param in_enriched the enrichment object
#' 
#' @export
#' @return data.table
enriched_to_fgsea = function(in_enriched)
{
  ranks = in_enriched@ranks
  results = data.table::as.data.table(in_enriched@statistics@statistic_data)
  rename_list = c(pval = "p", padj = "padjust", log2err = "log2err", ES = "es", NES = "nes", size = "size", leadingEdge = "leading_edge")
  results = dplyr::rename(results, dplyr::all_of(rename_list))
  results$pathway = in_enriched@statistics@annotation_id

  return(list(pathways = in_enriched@annotation@annotation_features,
             ranks = ranks,
             results = results))
}