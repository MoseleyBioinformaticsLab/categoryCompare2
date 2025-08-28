#' hypergeometric feature class
#'
#' class to hold features undergoing hypergeometric enrichment
#'
#' @slot significant the significant features
#' @slot universe all of the features measured
#' @slot annotation annotation object
#'
#' @export
setClass(
  "hypergeometric_features",
  slots = list(significant = "ANY", universe = "ANY", annotation = "annotation")
)

#' hypergeometric_features constructor
#'
#' function to easily construct a hypergeometric_features-class object.
#'
#' @param significant the features that are significant
#' @param universe the background or universe of features measured
#' @param annotation the annotation object
#'
#' @export
#' @return hypergeometric_features-class
hypergeometric_features = function(significant, universe, annotation) {
  new(
    "hypergeometric_features",
    significant = significant,
    universe = universe,
    annotation = annotation
  )
}

#' do hypergeometric enrichment
#'
#' @param hypergeometric_features a hypergeometric_features object
#' @param direction which direction to do the enrichment (over or under)
#' @param p_adjust how to correct the p-values (default is "BH")
#' @param min_features how many features should be annotated before testing it?
#'
#' @details The \code{min_features} argument here applies to the minumum number of features an annotation has
#'   from the universe of features supplied, \bold{not} the minumum number of features from the differential
#'   list.
#'   For more about the p-value adjustment, see \code{stats::p.adjust}
#' @export
#' @return enriched_result
#'
hypergeometric_feature_enrichment <- function(
  hypergeometric_features,
  direction = "over",
  p_adjust = "BH",
  min_features = 1
) {
  # cleanup the features and annotations (should be in separate function)
  hypergeometric_features@universe <- unique(hypergeometric_features@universe)

  tmp_annot_feature <- hypergeometric_features@annotation@annotation_features
  annotation_universe <- unique(unlist(tmp_annot_feature))

  hypergeometric_features@universe <- intersect(
    hypergeometric_features@universe,
    annotation_universe
  )
  tmp_annot_feature <- lapply(
    tmp_annot_feature,
    intersect,
    hypergeometric_features@universe
  )

  n_feature <- sapply(tmp_annot_feature, length)
  keep_annot <- n_feature >= min_features

  tmp_annot_feature <- tmp_annot_feature[keep_annot]

  hypergeometric_features@significant <- intersect(
    hypergeometric_features@significant,
    hypergeometric_features@universe
  )
  hypergeometric_features@annotation@annotation_features <- tmp_annot_feature

  # this probably needs its own function eventually
  if (length(hypergeometric_features@annotation@description) != 0) {
    hypergeometric_features@annotation@description <- hypergeometric_features@annotation@description[names(
      tmp_annot_feature
    )]
  }

  if (length(hypergeometric_features@annotation@links) != 0) {
    hypergeometric_features@annotation@links <- hypergeometric_features@annotation@links[names(
      tmp_annot_feature
    )]
  }

  # now get the counts annotated
  num_white_drawn <- sapply(
    hypergeometric_features@annotation@annotation_features,
    function(x) sum(hypergeometric_features@significant %in% x)
  )

  if (length(num_white_drawn) == 0) {
    num_white_drawn <- 0
  }

  num_white <- Biobase::listLen(
    hypergeometric_features@annotation@annotation_features
  )

  if (length(num_white) == 0) {
    num_white <- 0
  }

  num_black <- length(hypergeometric_features@universe) - num_white
  num_drawn <- length(hypergeometric_features@significant)

  hyper_stats <- hypergeometric_basic(
    num_white,
    num_black,
    num_drawn,
    num_white_drawn,
    direction
  )

  hyper_stats$counts <- num_white_drawn[names(hyper_stats$p)]
  hyper_stats$padjust <- p.adjust(hyper_stats$p, p_adjust)
  hyper_stats$significant <- rep(num_drawn, length(hyper_stats$padjust))

  out_stats <- new(
    "statistical_results",
    statistic_data = hyper_stats,
    annotation_id = names(
      hypergeometric_features@annotation@annotation_features
    ),
    method = "hypergeometric"
  )

  out_enrich <- new(
    "hypergeometric_result",
    features = hypergeometric_features@significant,
    universe = hypergeometric_features@universe,
    statistics = out_stats,
    annotation = hypergeometric_features@annotation
  )

  out_enrich
}

#' generate link text
#'
#' given a named vector of links, generate an actual html link formatted for output in html documents
#'
#' @param links the vector of links
#'
#' @export
#' @return character
generate_link <- function(links) {
  link_names <- names(links)

  if (max(nchar(link_names)) == 0) {
    link_names <- "link"
  }

  paste('<a href="', links, '">', link_names, '</a>', sep = "")
}

#' do hypergeometric test
#'
#' does a hypergeometric enrichment test
#'
#' @param num_white number of white balls in urn
#' @param num_black number of black balls in urn
#' @param num_drawn number of balls taken from urn
#' @param num_white_drawn number of white balls taken from urn
#' @param direction which direction is the test
#'
#' @export
#' @return list
hypergeometric_basic <- function(
  num_white,
  num_black,
  num_drawn,
  num_white_drawn,
  direction = "over"
) {
  n_2_1 <- num_white - num_white_drawn
  n_1_2 <- num_drawn - num_white_drawn
  n_2_2 <- num_black - n_1_2

  odds_ratio <- (num_white_drawn * n_2_2) / (n_1_2 * n_2_1)

  expected <- ((num_white_drawn + n_1_2) * (num_white_drawn + n_2_1)) /
    (num_white_drawn + n_1_2 + n_2_1 + n_2_2)

  p_values <- switch(
    direction,
    over = phyper(
      num_white_drawn - 1L,
      num_white,
      num_black,
      num_drawn,
      lower.tail = FALSE
    ),
    under = phyper(
      num_white_drawn,
      num_white,
      num_black,
      num_drawn,
      lower.tail = TRUE
    )
  )

  list(p = p_values, odds = odds_ratio, expected = expected)
}
