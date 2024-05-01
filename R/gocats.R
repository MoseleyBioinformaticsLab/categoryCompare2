#' gocats to annnotations
#' 
#' Transforms a gocats ancestors JSON list to a GO annotation object.
#' 
#' @param ancestors_file the ancestors.json file from gocats (required)
#' @param namespace_file the namespace.json file from gocats (optional)
#' @param annotation_type what annotations are we making? (gocatsGO by default)
#' @param feature_type what type of features are we using (assume Uniprot)
#' @param feature_translation a data.frame used to convert the feature IDs
#' 
#' @return annotation object
#' @export
gocats_to_annotation = function(ancestors_file = "ancestors.json",
                               namespace_file = "namespace.json",
                               annotation_type = "gocatsGO",
                               feature_type = "Uniprot",
                               feature_translation = NULL)
{
  stopifnot(file.exists(ancestors_file))
  
  ancestors = jsonlite::fromJSON(ancestors_file)
  
  if (!is.null(feature_translation)) {
    if (!inherits(feature_translation, "data.frame")) {
      stop("feature_translation must be a data.frame!")
    }
    if (!all(names(feature_translation) %in% c("from", "to"))) {
      stop("feature_translation must contain the columns 'from' and 'to'!")
    } else {
      match_names = intersect(feature_translation$from, names(ancestors))
      
      ancestors = ancestors[match_names]
      feature_translation = feature_translation[feature_translation$from %in% match_names, ]
      translations = feature_translation$to
      names(translations) = feature_translation$from
      translations = translations[match_names]
      names(ancestors) = translations
    }
    
  }
  
  go_2_gene = Biobase::reverseSplit(ancestors)
  go_2_gene = purrr::map(go_2_gene, unique)
  
  if (!is.null(namespace_file)) {
    if (!file.exists(namespace_file)) {
      warning(paste0(namespace_file, " does not exist. GO namespace will not be updated."))
      namespaces_short = character(0)
    } else {
      namespaces = jsonlite::fromJSON(namespace_file) |> unlist()
      namespaces_short = gsub("biological_process", "BP", namespaces)
      namespaces_short = gsub("molecular_function", "MF", namespaces_short)
      namespaces_short = gsub("cellular_component", "CC", namespaces_short)
    }
  } else {
    namespaces_short = character(0)
  }
  
  
  if (requireNamespace("GO.db", quietly = TRUE)) {
    descriptions = suppressMessages(AnnotationDbi::select(GO.db::GO.db, keys = names(go_2_gene), columns = "TERM", keytype = "GOID")$TERM)
    names(descriptions) = names(go_2_gene)
  } else {
    message("GO.db is not installed, no descriptions will be added to the GO terms.")
    descriptions = character(0)
  }
  
  if ((length(namespaces_short) > 0) && (length(descriptions) > 0)) {
    namespaces_short = namespaces_short[names(go_2_gene)]
    descriptions = descriptions[names(go_2_gene)]
    descriptions = paste0(namespaces_short, ":", descriptions)
    names(descriptions) = names(go_2_gene)
  }
  
  out_annotation = annotation(annotation_features = go_2_gene,
                              annotation_type = annotation_type,
                              description = descriptions,
                              feature_type = feature_type)
  return(out_annotation)
  
}
