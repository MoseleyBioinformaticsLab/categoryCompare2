#' binomial feature class
#'
#' class to hold features undergoing binomail statistical testing
#'
#' @slot positivefc the features with positive fold-changes
#' @slot negativefc the features with negative fold-changes
#' @slot annotation annotation object
#'
#' @export
setClass(
  "binomial_features",
  slots = list(
    positivefc = "ANY",
    negativefc = "ANY",
    annotation = "annotation"
  )
)

#' do binomial testing
#'
#' @param binomial_features a binomial_features object
#' @param p_expected the expected probability (default 0.5)
#' @param direction which direction to do the enrichment (two.sided, less, greater)
#' @param p_adjust how to correct the p-values (default is "BH")
#' @param conf_level the confidence level for the confidence interval (default is 0.95)
#' @param min_features a minimum number of features that are annotated to each annotation (default is 6)
#' @export
#' @return enriched_result
#'
binomial_feature_enrichment = function(
  binomial_features,
  p_expected = 0.5,
  direction = "two.sided",
  p_adjust = "BH",
  conf_level = 0.95,
  min_features = 6
) {
  # cleanup the features and annotations (should be in separate function)
  posfc = unique(binomial_features@positivefc)
  negfc = unique(binomial_features@negativefc)
  allfc = union(posfc, negfc)

  tmp_annot_feature = binomial_features@annotation@annotation_features
  annotation_universe = unique(unlist(tmp_annot_feature))
  annotation_universe = intersect(allfc, annotation_universe)

  posfc = intersect(posfc, annotation_universe)
  negfc = intersect(negfc, annotation_universe)

  binomial_features@positivefc = posfc
  binomial_features@negativefc = negfc
  tmp_annot_feature = lapply(tmp_annot_feature, intersect, annotation_universe)

  n_feature = sapply(tmp_annot_feature, length)
  keep_annot = n_feature >= min_features

  if (sum(keep_annot) == 0) {
    warning("No annotations had more than min_features, no tests ran!")
    return(NULL)
  }

  tmp_annot_feature = tmp_annot_feature[keep_annot]

  binomial_features@annotation@annotation_features = tmp_annot_feature

  # this probably needs its own function eventually
  if (length(binomial_features@annotation@description) != 0) {
    binomial_features@annotation@description = binomial_features@annotation@description[names(
      tmp_annot_feature
    )]
  }

  if (length(binomial_features@annotation@links) != 0) {
    binomial_features@annotation@links = binomial_features@annotation@links[names(
      tmp_annot_feature
    )]
  }

  # now get the number of positive and negative in each annotation
  feature_positive = lapply(
    binomial_features@annotation@annotation_features,
    function(x) {
      binomial_features@positivefc[binomial_features@positivefc %in% x]
    }
  )
  feature_negative = lapply(
    binomial_features@annotation@annotation_features,
    function(x) {
      binomial_features@negativefc[binomial_features@negativefc %in% x]
    }
  )
  num_positive = sapply(feature_positive, length)
  num_negative = sapply(feature_negative, length)

  binom_stats = binomial_basic(
    num_positive,
    num_positive + num_negative,
    p_expected,
    direction,
    conf_level = conf_level
  )

  binom_stats$padjust = p.adjust(binom_stats$p, p_adjust)
  binom_stats$num_positive = num_positive
  binom_stats$num_negative = num_negative

  direction = rep(0, length(num_positive))
  direction[num_positive > num_negative] = 1
  direction[num_positive < num_negative] = -1
  binom_stats$direction = direction

  out_stats = new(
    "statistical_results",
    statistic_data = binom_stats,
    annotation_id = names(binomial_features@annotation@annotation_features),
    method = "binomial"
  )

  out_enrich = new(
    "binomial_result",
    positivefc = posfc,
    negativefc = negfc,
    statistics = out_stats,
    annotation = binomial_features@annotation
  )

  out_enrich
}


#' do binomial test
#'
#' does a binomial test
#'
#' @param positive_cases number of positive instances
#' @param total_cases total number of cases observed
#' @param p_expected what is the expected probability
#' @param direction which direction is the test
#' @param conf_level confidence level for the confidence interval
#'
#' @export
#' @return list
binomial_basic = function(
  positive_cases,
  total_cases,
  p_expected = 0.5,
  direction = "two.sided",
  conf_level = 0.95
) {
  positive_cases = check_is_positive_integer(positive_cases, "positive_cases")
  total_cases = check_is_positive_integer(total_cases, "total_cases")

  if (length(positive_cases) != length(total_cases)) {
    stop("Number of values in positive_cases and total_cases does not match!")
  }

  if (
    (length(p_expected) != length(positive_cases)) && (length(p_expected) != 1)
  ) {
    stop(
      "p_expected should either match the number of positive_cases, or be a single number!"
    )
  }

  if (length(p_expected) == 1) {
    p_expected = rep(p_expected, length(positive_cases))
  }

  if ((min(p_expected) < 0) || (max(p_expected) > 1)) {
    stop("'p_expected' must be a value between 0 and 1")
  }

  less_binom = function(positive_cases, total_cases, p_expected) {
    stats::pbinom(positive_cases, total_cases, p_expected)
  }
  greater_binom = function(positive_cases, total_cases, p_expected) {
    stats::pbinom(
      positive_cases - 1,
      total_cases,
      p_expected,
      lower.tail = FALSE
    )
  }
  two_sided_binom = function(positive_cases, total_cases, p_expected) {
    two_sided_inner = function(positive_cases, total_cases, p_expected) {
      if ((p_expected == 0) && (positive_cases == 0)) {
        return(1)
      }
      if ((p_expected == 1) && (positive_cases == total_cases)) {
        return(0)
      }
      rel_err = 1 + 1e-07
      d_binom = stats::dbinom(positive_cases, total_cases, p_expected)
      m_values = total_cases * p_expected
      if (positive_cases == m_values) {
        return(1)
      } else if (positive_cases < m_values) {
        i_sequence = seq.int(from = ceiling(m_values), to = total_cases)
        y_values = sum(
          stats::dbinom(i_sequence, total_cases, p_expected) <=
            d_binom * rel_err
        )
        p_value = stats::pbinom(positive_cases, total_cases, p_expected) +
          stats::pbinom(
            total_cases - y_values,
            total_cases,
            p_expected,
            lower.tail = FALSE
          )
        return(p_value)
      } else {
        i_sequence = seq.int(from = 0, to = floor(m_values))
        y_values = sum(
          stats::dbinom(i_sequence, total_cases, p_expected) <=
            d_binom * rel_err
        )
        p_value = stats::pbinom(y_values - 1, total_cases, p_expected) +
          stats::pbinom(
            positive_cases - 1,
            total_cases,
            p_expected,
            lower.tail = FALSE
          )
        return(p_value)
      }
    }

    run_list = list(
      positive_cases = positive_cases,
      total_cases = total_cases,
      p_expected = p_expected
    )
    p_values = purrr::pmap_dbl(run_list, two_sided_inner)
    p_values
  }
  p_values = switch(
    direction,
    less = less_binom(positive_cases, total_cases, p_expected),
    greater = greater_binom(positive_cases, total_cases, p_expected),
    two.sided = two_sided_binom(positive_cases, total_cases, p_expected)
  )

  pvalues_lower = function(positive_cases, conf_level) {
    out_values = stats::qbeta(
      conf_level,
      positive_cases,
      total_cases - positive_cases + 1
    )
    out_values[positive_cases == 0] = 0
    return(out_values)
  }
  pvalues_upper = function(positive_cases, conf_level) {
    out_values = stats::qbeta(
      1 - conf_level,
      positive_cases + 1,
      total_cases - positive_cases
    )
    out_values[positive_cases == 1] = 1
    return(out_values)
  }
  pvalues_two_sided = function(positive_cases, conf_level) {
    conf_level = (1 - conf_level) / 2
    out_matrix = cbind(
      pvalues_lower(positive_cases, conf_level),
      pvalues_upper(positive_cases, conf_level)
    )
    return(out_matrix)
  }
  conf_interval = switch(
    direction,
    less = cbind(
      rep(0, length(positive_cases)),
      pvalues_upper(positive_cases, 1 - conf_level)
    ),
    greater = cbind(
      pvalues_lower(positive_cases, 1 - conf_level),
      rep(1, length(positive_cases))
    ),
    two.sided = pvalues_two_sided(positive_cases, conf_level)
  )

  return(list(
    statistic = positive_cases,
    parameter = total_cases,
    p = p_values,
    conf_lower = conf_interval[, 1],
    conf_upper = conf_interval[, 2],
    estimate = positive_cases / total_cases,
    null_value = p_expected,
    alternative = direction
  ))
}


check_is_positive_integer = function(values, variable_name) {
  if (any(values < 0)) {
    stop(paste0(
      variable_name,
      " contains negative values! This is not allowed for this test."
    ))
  }

  if (any(is.na(values))) {
    stop(paste0(
      variable_name,
      " contains NA values! This is not allowed for this test."
    ))
  }

  if (any(is.infinite(values))) {
    stop(paste0(
      variable_name,
      " contains Inf values! This is not allowed for this test."
    ))
  }

  values_round = round(values)
  round_diff = abs(values - values_round)
  if (max(round_diff) > 1e-07) {
    stop(paste0(
      variable_name,
      " contains values that are non-integer! This is not allowed for this test."
    ))
  }
  values_round
}
