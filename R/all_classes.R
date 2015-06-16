#' annotation class
#' 
#' This class holds an annotation object that defines how annotations relate to
#' features, as well as various pieces about each annotation
#' 
#' These objects may be created by hand, or may result from specific functions
#' to create them. Most notably, this package provides functions for creating
#' a Gene Ontology annotation.
#' 
#' @slot annotation_features list of annotation to feature relationships
#' @slot description character vector providing descriptive text about the annotation
#' @slot counts numeric vector of how many features are in each annotation
#' @slot links character vector defining html links for each annotation (may be empty)
#' @slot type a one word short description of the "type" of annotation
#' 
#' @export
annotation <- setClass("annotation",
                       slots = list(annotation_features = "list",
                                    description = "character",
                                    counts = "numeric",
                                    links = "character",
                                    type = "character"))

#' show annotation
#' @exportMethod show
setMethod("show", signature = list(object = "annotation"),
          function(object){
            n_gene <- length(unique(unlist(object@annotation_features)))
            n_annot <- length(object@annotation_features)
            cat(object@type, "Annotation Object\n", sep = " ")
            cat("with ", n_annot, " annotations and ", n_gene, " genes\n", sep = "")
          })

#' annotation constructor
#' 
#' Does sensical checks when creating an \linkS4class{annotation} object.
#' 
#' See the \linkS4class{annotation}, each slot is a parameter.
#' 
#' @param annotation_features list of annotation to feature relationships
#' @param description character vector providing descriptive text about the annotation
#' @param counts numeric vector of how many features are in each annotation
#' @param links character vector defining html links for each annotation (may be empty)
#' @param type a one word short description of the "type" of annotation
#' 
#' @rdname annotation_constructor
#' @export
annotation <- function(annotation_features, type = NULL, description = character(0), links = character(0)){
  if (is.null(type)){
    type <- "generic_annotation"
  }
  
  annotation_names <- names(annotation_features)
  n_annot <- length(annotation_names)
  
  if (length(description) > 0){
    desc_names <- names(description)
    if (sum(annotation_names %in% desc_names) == n_annot){
      description <- description[annotation_names]
    } else {
      stop("description names don't match annotation names!", call. = FALSE)
    }
  }
  
  if (length(links) > 0){
    link_names <- names(links)
    if (sum(annotation_names %in% link_names) == n_annot){
      links <- links[annotation_names]
    } else {
      stop("description names don't match annotation names!", call. = FALSE)
    }
  }
  
  counts <- vapply(annotation_features, length, numeric(1))
  
  new("annotation",
      type = type,
      description = description,
      links = links,
      counts = counts)
  
}

#' statistical results class
#' 
#' This class holds the part of an enrichment that is the statistical results.
#' It has two pieces, a \code{list} of \code{statistics} that is a named list
#' with the actual numerical results of applying the statistics. The other piece
#' is the \code{annotation_id} vector defining which entry in each vector of the
#' \code{statistics} is.
#' 
#' @slot statistic_data list of numerical statistics
#' @slot annotation_id vector of ids
#' 
#' @export
statistical_results <- setClass("statistical_results",
                                slots = list(statistic_data = "list",
                                             annotation_id = "ANY"))

#' the enriched results class
#' 
#' @slot features the "features" of interest, a vector of class ANY
#' @slot universe all of the "features" in the background
#' @slot annotation list giving the annotation to feature relationship
#' @slot statistics a \linkS4class{statistical_results} object
#' 
#' @export
enriched_result <- setClass("enriched_result",
                            slots = list(features = "ANY",
                                         universe = "ANY",
                                         annotation = "annotation",
                                         statistics = "statistical_results"))

#' creates enriched result
#' 
#' given all the slots for an \linkS4class{enriched_result}, checks that all
#' the data is self-consistent, and creates the \code{enriched_result} object.
#' 
#' @param features the features that were differentially expressed (see details)
#' @param universe all of the features that were measured
#' @param annotation an \linkS4class{annotation} object
#' @param statistics a \linkS4class{statistical_results} object
#' 
#' @export
#' @return enriched_result
#' @rdname enriched_result_constructor
enriched_result <- function(features, universe, annotation, statistics){
  stopifnot(class(statistics) == "statistical_results")
  stopifnot(class(annotation) == "annotation")
  
  annotation_names <- names(annotation@annotation_features)
  stat_names <- statistics@annotation_id
  
  annotation_stat <- intersect(annotation_names, stat_names)
  
  if (sum(stat_names %in% annotation_stat) != length(stat_names)){
    stop("There are missing annotations!", call. = TRUE)
  }
  
  annotation@annotation_features <- annotation@annotation_features[stat_names]
  if (length(annotation@description) != 0){
    annotation@description <- annotation@description[stat_names]
  }
  
  if (length(annotation@counts) != 0){
    annotation@counts <- annotation@counts[stat_names]
  }
  
  if (length(annotation@links) != 0){
    annotation@links <- annotation@links[stat_names]
  }
  
  return(new("enriched_result",
             features = features,
             universe = universe,
             annotation = annotation,
             statistics = statistics))
}


#' significant annotations
#' 
#' The \code{significant_annotations} class holds which annotations from which
#' enrichment were both \strong{measured} and \strong{significant}. Each of these
#' slots is a \emph{logical matrix} with rows named by \emph{annotation_id} and 
#' columns named by the names of the \linkS4class{enriched_result} that was combined.
#' 
#' @slot significant logical matrix
#' @slot measured logical matrix
#' @slot sig_calls character representations of calls used to filter the data
#' 
#' @export
significant_annotations <- setClass("significant_annotations",
                                    slots = list(significant = "matrix",
                                                 measured = "matrix",
                                                 sig_calls = "character"))

#' create significant annotations
#' 
#' Makes a new significant_annotation while checking that everything is valid.
#' 
#' @param significant logical matrix of annotations (rows) and experiments (columns)
#' @param measured logical matrix of annotations (rows) and experiments (columns)
#' @param sig_calls character vector of deparsed calls that resulted in signficant and
#' measured
#' 
#' @export
significant_annotations <- function(significant, measured, sig_calls = NULL){
  stopifnot(nrow(significant) != nrow(measured))
  stopifnot(ncol(significant) != ncol(measured))
  
  new("signficant_annotations",
      significant = significant,
      measured = measured,
      sig_calls = sig_calls)
}

#' combined statistics
#' 
#' holds the results of extracting a bunch of statistics from a \linkS4class{combined_enrichment}
#' into one entity. This is useful because we want to enable multiple data representations and
#' simple filtering on the actual \code{data.frame} of statistics, and this provides flexibility
#' to enable that.
#' 
#' @slot statistic_data a \code{data.frame} of all of the statistics from all of the enrichments
#' @slot significant a \linkS4class{significant_annotations} object, that may be empty
#' @slot which_enrichment a \code{vector} giving which enrichment each column of the statistics came from
#' @slot which_statistic a \code{vector} providing which statistic each column contains
#' 
#' @export
combined_statistics <- setClass("combined_statistics",
                                contains = "statistical_results",
                                slots = list(statistic_data = "data.frame",
                                             annotation_id = "character",
                                             significant = "significant_annotations",
                                             which_enrichment = "character",
                                             which_statistic = "character"))

#' combined enrichments
#' 
#' The \code{combined_enrichment} class holds the results of combining several 
#' \linkS4class{enriched_result}s together, which includes the original 
#' \linkS4class{enriched_result}s, as well as the \code{annotation_graph}
#' and combined \linkS4class{annotation} objects.
#' 
#' @slot enriched list of enriched objects
#' @slot enriched_type character describing the enrichment annotation
#' @slot annotation \linkS4class{annotation} where the annotation_features
#' have been combined across the \linkS4class{enriched_results}
#' 
#' @export
combined_enrichment <- setClass("combined_enrichment",
                                slots = list(enriched = "list",
                                             annotation = "annotation",
                                             statistics = "combined_statistics"))

#' cc_graph
#' 
#' A \code{cc_graph} class is a \code{graphNEL} with the added slot of
#' \code{significant}, a matrix of rows (nodes / annotations) and whether
#' they were found to be significant in a given enrichment (columns). This
#' matrix is used for classifying the annotations into different groups, and
#' generating either pie-charts or coloring the nodes in a visualization.
#' 
#' @slot significant numeric matrix of ones and zeros
#' 
#' @export
cc_graph <- setClass("cc_graph",
                     contains = "graphNEL",
                     slots = list(significant = "matrix"))

#' cc_graph constructor
#' 
#' constructs a \emph{cc_graph} given a \linkS4class{graphNEL} and a \emph{significant} matrix.
#' 
#' @param graph the \linkS4class{graphNEL}
#' @param significant a matrix indicating which nodes are significant in which experiment
#' 
#' @export
cc_graph <- function(graph, significant){
  out_graph <- as(graph, "cc_graph")
  out_graph@significant <- significant
  out_graph
}

#' node_assign
#' 
#' The \code{node_assign} class holds the unique annotation combinations and the
#' assignment of the nodes to those combinations for use in visualization.
#' 
#' @slot groups the unique groups, as a logical matrix
#' @slot assignments named character vector providing association with groups
#' @slot colors named character vector of hex colors for groups or experiments
#' @slot color_type whether doing group or experiment based colors
#' @slot pie_locs if doing experiment colors, then pie graphs were generated here
#' 
#' @export
node_assign <- setClass("node_assign",
                        slots = list(groups = "matrix",
                                     assignments = "character",
                                     colors = "character",
                                     color_type = "character",
                                     pie_locs = "character"))
