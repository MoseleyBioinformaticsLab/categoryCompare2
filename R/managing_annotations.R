#' orgdb annotations
#' 
#' Generate an annotation object for genes based on an "org.*.db" object, and pulling
#' information from it.
#' 
#' @param orgdb the name of the org.*.db object
#' @param feature_type which type of IDs to map (see details)
#' @param annotation_type the type of annotation to grab (see details)
#' 
#' @md
#' 
#' @details This function generates a `categoryCompare2` annotation object
#' from a Bioconductor "org.*.db" object. Even though different gene identifiers can
#' be used, almost all of the mappings are via ENTREZID.
#' 
#' The set of feature or gene keys that can be used to create the annotations include:
#' * ENTREZID: ENTREZ gene ids
#' * ACCNUM: genbank accession numbers
#' * SYMBOL: gene symbols, eg ABCA1
#' * GENENAME: gene names, eg "ATP binding cassette subfamily A member 1"
#' * ENSEMBL: the ensembl gene ids (all start with ENSG...)
#' * ENSEMBLPROT: ensembl protein ids (ENSP...)
#' * ENSEMBLTRANS: ensemlb transcript ids (ENST...)
#' * REFSEQ: reference sequence IDs, NM, NP, NR, XP, etc
#' * UNIGENE: gene ids from UNIPROT eg Hs.88556
#' * UNIPROT: protein ids from UNIPROT eg P80404
#' 
#' The set of annotations that can be mapped to features include:
#' * GO: annotations from gene ontology
#' * PATH: KEGG Pathway identifiers (not updated since 2011!)
#' * CHRLOC: location on the chromosome
#' * OMIM: mendelian inheritance in man identifiers
#' * PMID: pubmed identifiers
#' * PROSITE
#' * PFAM: protein family identifiers
#' * IPI: protein-protein interactions
#' 
#' For GO annotations, it is also possible to pass `GO` to use all 3 sub-ontologies simultaneously,
#' or any combination of `BP`, `MF`, and `CC`.
#' 
#' @export
#' @return annotation object
#' 
get_db_annotation <- function(orgdb = "org.Hs.eg.db", features = NULL, feature_type = "ENTREZID",
                              annotation_type = "GO"){
  go_types <- c("GO", "BP", "MF", "CC")
  go_sub <- c("BP", "MF", "CC")
  
  if (!require(orgdb, character.only = TRUE, quietly = TRUE)) {
    stop("The package ", orgdb, " is not installed/available. Try installing it with biocLite('", orgdb, "')")
  } else {
    annotation_src <- eval(parse(text = orgdb))
    annotation_columns <- AnnotationDbi::columns(annotation_src)
    annotation_keytypes <- AnnotationDbi::keytypes(annotation_src)
    
    if (!(annotation_type %in% c(go_types, annotation_columns))) {
      stop("Unknown annotation type!")
    }
    
    if (!(feature_type %in% annotation_keytypes)) {
      stop("Unknown feature_type!")
    }
    
    if (is.null(features)) {
      features <- AnnotationDbi::keys(annotation_src, feature_type)
    }
    
    if (annotation_type %in% go_types) {
      suppressMessages(require("GO.db", character.only = TRUE))
      feature_ann_map <- suppressMessages(AnnotationDbi::select(annotation_src, keys = features,
                                                             keytype = feature_type,
                                                             columns = "GOALL"))
      
      if (annotation_type %in% go_sub) {
        feature_ann_map <- feature_ann_map[feature_ann_map$ONTOLOGYALL %in% annotation_type, ]
      }
      
      
      ann_feature_list <- split(feature_ann_map[[feature_type]], feature_ann_map[["GOALL"]])
      ann_feature_list <- lapply(ann_feature_list, unique)
      ann_description <- suppressMessages(AnnotationDbi::select(GO.db, keys = names(ann_feature_list), columns = "TERM", keytype = "GOID")$TERM)
      names(ann_description) <- names(ann_feature_list)
      
      if (annotation_type %in% "GO") {
        go_ontology_map <- unique(feature_ann_map[, c("GOALL", "ONTOLOGYALL")])
        go_ontology <- go_ontology_map$ONTOLOGYALL
        names(go_ontology) <- go_ontology_map$GOALL
        go_ontology <- go_ontology[names(ann_description)]
        ann_description <- paste0(go_ontology, ":", ann_description)
        names(ann_description) <- names(go_ontology)
      }
      
      
      annotation_obj <- categoryCompare2::annotation(annotation_features = ann_feature_list,
                                                     description = ann_description,
                                                     annotation_type = annotation_type,
                                                     feature_type = feature_type)
    } else {
      feature_ann_map <- suppressMessages(AnnotationDbi::select(annotation_src, keys = features,
                                                                keytype = feature_type,
                                                                columns = annotation_type))
      ann_feature_list <- split(feature_ann_map[[feature_type]], feature_ann_map[[annotation_type]])
      ann_feature_list <- lapply(ann_feature_list, unique)
      
      annotation_obj <- categoryCompare2::annotation(annotation_features = ann_feature_list,
                                                     annotation_type = annotation_type,
                                                     feature_type = feature_type)
    }
    
  }
  annotation_obj
}


#' annotation to json
#' 
#' Given a `categoryCompare2` annotation object, generate a JSON representation
#' that can be used with the command line executable
#' 
#' @param annotation_obj the annotation object
#' @param json_file the file to save it to
#' 
#' @return the json string (invisibly)
#' @export
annotation_2_json <- function(annotation_obj, json_file = NULL){
  obj_list <- purrr::map(slotNames(annotation_obj), function(x){
    tmp_data <- slot(annotation_obj, x)
    
    if (length(tmp_data) != 0) {
      if (x %in% c("description", "links")) {
        out_data <- as.list(tmp_data)
      } else {
        out_data <- tmp_data
      }
    } else {
      out_data <- NULL
    }
    out_data
  })
  names(obj_list) <- slotNames(annotation_obj)
  
  obj_list <- obj_list[!purrr::map_lgl(obj_list, is.null)]
  
  obj_json <- jsonlite::toJSON(obj_list, pretty = TRUE, auto_unbox = TRUE)
  if (!is.null(json_file)) {
    cat(obj_json, file = json_file, sep = "\n")
  } 
  invisible(obj_json)
}

replace_null <- function(x){
  if (is.null(x)) {
    NA
  } else {
    x
  }
}

#' json to annotation
#' 
#' Given a JSON based annotation object, read it in and create the `annotation`
#' for actually doing enrichment. 
#' 
#' @param json_file the json annotation file
#' 
#' @return annotation object
#' @export
json_2_annotation <- function(json_file){
  stopifnot(file.exists(json_file))
  annotation_list <- jsonlite::fromJSON(json_file, simplifyVector = TRUE)
  #print(names(annotation_list))
  
  if (is.null(names(annotation_list))) {
    annotation_list <- annotation_list[[1]]
  }
  
  if (!is.null(annotation_list$description)) {
    annotation_list$description <- purrr::map(annotation_list$description, replace_null)
    annotation_list$description <- unlist(annotation_list$description, use.names = TRUE)
  } else {
    annotation_list$description <- character(0)
  }
  if (!is.null(annotation_list$links)) {
    annotation_list$links <- unlist(annotation_list$links)
  } else {
    annotation_list$links <- character(0)
  }

  #print(names(annotation_list))
  #print(annotation_list$description)
  annotation(annotation_features = annotation_list$annotation_features,
             annotation_type = annotation_list$annotation_type,
             description = annotation_list$description,
             links = annotation_list$links,
             feature_type = annotation_list$feature_type)
}

#' annotation reversal
#' 
#' Given a JSON file of features to annotations, reverse to turn it into
#' annotations to features, and optionally add some meta-information about them.
#' 
#' @param json_file the json file to use
#' @param out_file the json file to write out to
#' @param feature_type the type of features
#' @param annotation_type the type of annotations
#' 
#' @importFrom jsonlite fromJSON toJSON
#' @export
#' @return the json object, invisibly
#' 
json_annotation_reversal <- function(json_file, out_file = "annotations.json",
                                feature_type = NULL, annotation_type = NULL){
  stopifnot(file.exists(json_file))
  
  in_annotation <- jsonlite::fromJSON(json_file, simplifyVector = FALSE, flatten = TRUE)
  if (length(in_annotation) == 1) {
    in_annotation <- in_anntation[[1]]
  }
  
  rev_annotation <- Biobase::reverseSplit(in_annotation)
  rev_annotation <- purrr::map(rev_annotation, unique)
  
  out_annotation <- annotation(annotation_features = rev_annotation,
                               description = character(0),
                               links = character(0),
                               annotation_type = annotation_type,
                               feature_type = feature_type)
  
  out_json <- annotation_2_json(out_annotation, out_file)
  out_json
}
