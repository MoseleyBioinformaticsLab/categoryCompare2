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
#' @slot annotation_type a one word short description of the "type" of annotation
#' @slot feature_type a one word short description of the "type" of features
#' 
#' @export
#' @rdname annotation
annotation <- setClass("annotation",
                       slots = list(annotation_features = "list",
                                    description = "character",
                                    counts = "numeric",
                                    links = "character",
                                    annotation_type = "character",
                                    feature_type = "character"))

#' show annotation
#' 
#' @param object the annotation object
#' 
#' @exportMethod show
#' @rdname annotation
setMethod("show", signature = list(object = "annotation"),
          function(object){
            n_gene <- length(unique(unlist(object@annotation_features)))
            n_annot <- length(object@annotation_features)
            cat("      Annotation Type:", object@annotation_type, "\n")
            cat("         Feature Type:", object@feature_type, "\n")
            cat("Number of Annotations:", n_annot, "\n")
            cat("   Number of Features:", n_gene, "\n")
          })

#' annotation constructor
#' 
#' Does sensical checks when creating an \code{annotation} object.
#' 
#' See the \code{annotation}, each slot is a parameter.
#' 
#' @param annotation_features list of annotation to feature relationships
#' @param annotation_type a simple one word description of the annotations
#' @param description character vector providing descriptive text about the annotation
#' @param links character vector defining html links for each annotation (may be empty)
#' @param feature_type one word description of the feature type
#' 
#' @rdname annotation
#' @export
annotation <- function(annotation_features, annotation_type = NULL, description = character(0), links = character(0),
                       feature_type = NULL){
  if (is.null(annotation_type)){
    annotation_type <- "GENERIC"
  }
  if (is.null(feature_type)) {
    feature_type <- "UNKNOWN"
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
  
  # print(class(annotation_features))
  # print(class(annotation_type))
  # print(class(description))
  # print(class(links))
  # print(class(feature_type))
  # print(class(counts))
  
  new("annotation",
      annotation_features = annotation_features,
      annotation_type = annotation_type,
      feature_type = feature_type,
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
#' @slot method how the statistics were calculated
#' 
#' @export
statistical_results <- setClass("statistical_results",
                                slots = list(statistic_data = "list",
                                             annotation_id = "ANY",
                                             method = "character"))

#' the enriched results class
#' 
#' @slot features the "features" of interest, a vector of class ANY
#' @slot universe all of the "features" in the background
#' @slot annotation list giving the annotation to feature relationship
#' @slot statistics a \code{\link{statistical_results}} object
#' 
#' @export
enriched_result <- setClass("enriched_result",
                            slots = list(features = "ANY",
                                         universe = "ANY",
                                         ranks = "ANY",
                                         annotation = "annotation",
                                         statistics = "statistical_results"))

#' show enriched_result
#' 
#' @param object the enriched_result object to show
#' 
#' @exportMethod show
setMethod("show", signature = list(object = "enriched_result"),
          function(object){
            stat_method <- object@statistics@method
            n_feature <- length(unique(object@features))
            n_universe <- length(unique(object@universe))
            
            enrich_type <- object@annotation@annotation_type
            n_annot <- length(object@annotation@annotation_features)
            
            cat("   Enrichment Method:", stat_method, "\n")
            cat("     Annotation Type:", enrich_type, "\n")
            cat("Significant Features:", n_feature, "\n")
            cat("       Universe Size:", n_universe, "\n")
          })

#' creates enriched result
#' 
#' given all the slots for an \code{\link{enriched_result}}, checks that all
#' the data is self-consistent, and creates the \code{\link{enriched_result}} object.
#' 
#' @param features the features that were differentially expressed (see details)
#' @param universe all of the features that were measured
#' @param annotation an \code{\link{annotation}} object
#' @param statistics a \code{\link{statistical_results}} object
#' 
#' @export
#' @return enriched_result
#' @rdname enriched_result
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


#' the binomial results class
#' 
#' @slot positivefc the positive log-fold-changed genes, a vector of class ANY
#' @slot negativefc the negative log-fold-changed genes
#' @slot annotation list giving the annotation to feature relationship
#' @slot statistics a \code{\link{statistical_results}} object
#' 
#' @export
binomial_result <- setClass("binomial_result",
                            slots = list(positivefc = "ANY",
                                         negativefc = "ANY",
                                         annotation = "annotation",
                                         statistics = "statistical_results"))

#' show binomial_result
#' 
#' @param object the binomial_result object to show
#' 
#' @exportMethod show
setMethod("show", signature = list(object = "binomial_result"),
          function(object){
            stat_method <- object@statistics@method
            
            enrich_type <- object@annotation@annotation_type
            n_annot <- length(object@annotation@annotation_features)
            
            cat("   Enrichment Method:", stat_method, "\n")
            cat("     Annotation Type:", enrich_type, "\n")
          })

#' significant annotations
#' 
#' The \code{significant_annotations} class holds which annotations from which
#' enrichment were both \strong{measured} and \strong{significant}. Each of these
#' slots is a \emph{logical matrix} with rows named by \emph{annotation_id} and 
#' columns named by the names of the \code{\link{enriched_result}} that was combined.
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

#' show signficant_annotations
#' 
#' @param object the significant annotations object to show
#' 
#' @exportMethod show
setMethod("show", signature = list(object = "significant_annotations"),
          function(object){
            sig_calls <- object@sig_calls
            if (length(sig_calls) == 0){
              sig_calls <- "none"
            }
            out_calls <- paste("  ", sig_calls, sep = "")
            cat("Signficance Cutoffs:\n")
            cat(out_calls, sep = "\n")
            cat("\nCounts:\n")
            node_assign <- annotation_combinations(object)
            print(node_assign)
          })

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
#' holds the results of extracting a bunch of statistics from a \code{\link{combined_enrichment}}
#' into one entity. This is useful because we want to enable multiple data representations and
#' simple filtering on the actual \code{data.frame} of statistics, and this provides flexibility
#' to enable that.
#' 
#' @slot statistic_data a \code{data.frame} of all of the statistics from all of the enrichments
#' @slot significant a \code{\link{significant_annotations}} object, that may be empty
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

#' show combined_statistics
#' 
#' @param object \code{\link{combined_statistics}}
#' 
#' @exportMethod show
setMethod("show", signature = list(object = "combined_statistics"),
          function(object){
            print(object@significant)
          })

#' combined enrichments
#' 
#' The \code{combined_enrichment} class holds the results of combining several 
#' \code{\link{enriched_result}}s together, which includes the original 
#' \code{\link{enriched_result}}s, as well as the \code{\link{cc_graph}}
#' and combined \code{\link{annotation}} objects.
#' 
#' @slot enriched list of enriched objects
#' @slot annotation \code{\link{annotation}} where the annotation_features
#' have been combined across the \code{\link{enriched_result}}
#' @slot statistics \code{\link{combined_statistics}} of both
#' 
#' @export
#' @rdname combined_enrichment
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
#' @importFrom graph graphNEL
#' @rdname cc_graph
cc_graph <- setClass("cc_graph",
                     contains = "graphNEL",
                     slots = list(significant = "matrix"))

#' show cc_graph enrichment
#' 
#' @param object the cc_graph to show
#' 
#' @exportMethod show
#' @rdname cc_graph
setMethod("show", signature = list(object = "cc_graph"),
          function(object) {
            numNodes <- graph::numNodes(object)
            numEdge <- graph::numEdges(object)
            node_assign <- annotation_combinations(object)
            cat("A cc_graph with\n")
            cat("Number of Nodes =", numNodes, "\n")
            cat("Number of Edges =", numEdge, "\n")
            print(node_assign)
          })

#' cc_graph constructor
#' 
#' constructs a \emph{cc_graph} given a \code{\link[graph]{graphNEL}} and a \emph{significant} matrix.
#' 
#' @param graph the \code{\link[graph]{graphNEL}}
#' @param significant a matrix indicating which nodes are significant in which experiment
#' 
#' @export
#' @rdname cc_graph
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
#' @slot description named character vector providing a description to group
#' @slot colors named character vector of hex colors for groups or experiments
#' @slot color_type whether doing group or experiment based colors
#' @slot pie_locs if doing experiment colors, then pie graphs were generated here
#' 
#' @export
node_assign <- setClass("node_assign",
                        slots = list(groups = "matrix",
                                     assignments = "character",
                                     description = "character",
                                     colors = "character",
                                     color_type = "character",
                                     pie_locs = "ANY"))

#' show node_assign
#' 
#' @param object the node_assign to see
#' 
#' @exportMethod show
setMethod("show", signature = list(object = "node_assign"),
          function(object){
            group_matrix <- object@groups
            assignments <- object@assignments
            
            group_names <- rownames(group_matrix)
            numerical_group <- apply(group_matrix, 2, as.numeric)
            
            # we do this in case we have created a vector instead of a matrix from the 
            # apply step above
            numerical_group <- matrix(numerical_group, nrow = nrow(group_matrix), ncol = ncol(group_matrix))
            rownames(numerical_group) <- group_names
            colnames(numerical_group) <- colnames(group_matrix)
            
            group_counts <- vapply(group_names, function(x){sum(assignments %in% x)}, numeric(1))
            numerical_group <- cbind(numerical_group, counts=group_counts)
            print(numerical_group)
          })
