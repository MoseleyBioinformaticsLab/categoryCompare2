# runs the tests for getting significant results from combined objects

test_that("list indexing works properly", {
  c1 <- list(
    pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
    counts = c(a1 = 5, a2 = 10, a3 = 1),
    odds = c(a1 = 20, a2 = 100, a3 = 0)
  )

  pvalue1 <- c(TRUE, FALSE, TRUE) # pvalues < 0.05
  pvalue2 <- c(FALSE, FALSE, TRUE) # pvalues < 0.001

  counts1 <- c(TRUE, TRUE, FALSE) # counts > 1
  counts2 <- c(TRUE, FALSE, TRUE) # counts < 10

  odds1 <- c(TRUE, TRUE, FALSE) # odds > 0

  pvalue_counts1 <- c(TRUE, FALSE, TRUE) # pvalues < 0.5, counts >= 1

  # rlang::enquos doesn't work at the top level, so we need a little
  # wrapper function to be able to test this properly
  test_fun = function(list_to_query, ...) {
    queries = rlang::enquos(...)
    multi_query_list(list_to_query, queries)
  }

  expect_equal(
    pvalue1,
    test_fun(c1, pvalues < 0.05)
  )

  expect_equal(pvalue2, test_fun(c1, pvalues < 0.001))
  expect_equal(counts1, test_fun(c1, counts > 1))
  expect_equal(counts2, test_fun(c1, counts < 10))
  expect_equal(odds1, test_fun(c1, odds > 0))
  expect_equal(pvalue_counts1, test_fun(c1, pvalues < 0.5, counts >= 1))
  expect_error(test_fun(c1, asdfghjkl))
})


test_that("get correct significant annotations", {
  test_stat <- new(
    "statistical_results",
    annotation_id = c("a1", "a2", "a3"),
    statistic_data = list(
      pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
      counts = c(a1 = 5, a2 = 10, a3 = 1),
      odds = c(a1 = 20, a2 = 100, a3 = 0)
    )
  )

  expect_equal(
    c("a1", "a3"),
    get_significant_annotations(test_stat, pvalues < 0.05)
  )
  expect_equal(c("a1", "a2"), get_significant_annotations(test_stat, odds > 10))
  expect_equal(
    c("a1", "a3"),
    get_significant_annotations(test_stat, pvalues < 0.05, counts >= 1)
  )
})

# here we test what we get back from a combined_enrichment object

test_that("get correct sig annotations from combined_enrichment", {
  stat1 <- new(
    "statistical_results",
    annotation_id = c("a1", "a2", "a3"),
    statistic_data = list(
      pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
      counts = c(a1 = 5, a2 = 10, a3 = 1),
      odds = c(a1 = 20, a2 = 100, a3 = 0)
    )
  )
  stat2 <- new(
    "statistical_results",
    annotation_id = c("a1", "a2", "a4"),
    statistic_data = list(
      pvalues = c(a1 = 0.01, a2 = 0.03, a4 = 0.0001),
      counts = c(a1 = 5, a2 = 10, a4 = 1),
      odds = c(a1 = 20, a2 = 100, a4 = 0)
    )
  )

  en1 <- new(
    "hypergeometric_result",
    features = letters,
    universe = letters,
    annotation = new("annotation"),
    statistics = stat1
  )

  en2 <- new(
    "hypergeometric_result",
    features = letters,
    universe = letters,
    annotation = new("annotation"),
    statistics = stat2
  )

  test_combined <- new(
    "combined_enrichment",
    enriched = list(en1 = en1, en2 = en2),
    annotation = new("annotation")
  )

  combined_stats <- extract_statistics(test_combined)
  test_combined@statistics <- combined_stats

  # this uses pvalues < 0.05
  meas_matrix <- matrix(
    c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE),
    nrow = 4,
    ncol = 2
  )
  rownames(meas_matrix) <- c("a1", "a2", "a3", "a4")
  colnames(meas_matrix) <- c("en1", "en2")
  sig_matrix <- meas_matrix
  sig_matrix["a2", ] <- c(FALSE, TRUE)

  sig_call = "~pvalues < 0.05"
  names(sig_call) = ""
  expected_significant_annotations <- new(
    "significant_annotations",
    significant = sig_matrix,
    measured = meas_matrix,
    sig_calls = sig_call
  )

  expected_sig_combined <- test_combined
  expected_sig_combined@statistics@significant <- expected_significant_annotations

  expect_equal(
    expected_sig_combined,
    get_significant_annotations(test_combined, pvalues < 0.05)
  )
})


test_that("get correct sig annotations from combined_enrichment with single case", {
  stat1 <- new(
    "statistical_results",
    annotation_id = c("a1", "a2", "a3"),
    statistic_data = list(
      pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
      counts = c(a1 = 5, a2 = 10, a3 = 1),
      odds = c(a1 = 20, a2 = 100, a3 = 0)
    )
  )
  stat2 <- new(
    "statistical_results",
    annotation_id = c("a1", "a2", "a4"),
    statistic_data = list(
      pvalues = c(a1 = 0.01, a2 = 0.03, a4 = 0.0001),
      counts = c(a1 = 5, a2 = 10, a4 = 1),
      odds = c(a1 = 20, a2 = 100, a4 = 0)
    )
  )

  en1 <- new(
    "hypergeometric_result",
    features = letters,
    universe = letters,
    annotation = new("annotation"),
    statistics = stat1
  )

  en2 <- new(
    "hypergeometric_result",
    features = letters,
    universe = letters,
    annotation = new("annotation"),
    statistics = stat2
  )
  single_combined = new(
    "combined_enrichment",
    enriched = list(en1 = en1),
    annotation = new("annotation")
  )

  combined_stats <- extract_statistics(single_combined)
  single_combined@statistics <- combined_stats

  # this uses pvalues < 0.05
  meas_matrix <- matrix(
    c(TRUE, TRUE, TRUE),
    nrow = 3,
    ncol = 1
  )
  rownames(meas_matrix) <- c("a1", "a2", "a3")
  colnames(meas_matrix) <- c("en1")
  sig_matrix <- meas_matrix
  sig_matrix["a2", ] <- c(FALSE)

  sig_call = "~pvalues < 0.05"
  names(sig_call) = ""
  expected_significant_annotations <- new(
    "significant_annotations",
    significant = sig_matrix,
    measured = meas_matrix,
    sig_calls = sig_call
  )

  expected_single_combined <- single_combined
  expected_single_combined@statistics@significant <- expected_significant_annotations

  expect_equal(
    expected_single_combined,
    get_significant_annotations(single_combined, pvalues < 0.05)
  )
})
