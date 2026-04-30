test_that("various binomial flavors work", {
  # setup -----------
  h_10_neg = enrichment_data$table10 %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::pull(Entrez) %>%
    unique()

  h_10_pos = enrichment_data$table10 %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::pull(Entrez) %>%
    unique()

  h_48_neg = enrichment_data$table48 |>
    dplyr::filter(logFC < 0) |>
    dplyr::pull(Entrez) |>
    unique()

  h_48_pos = enrichment_data$table48 |>
    dplyr::filter(logFC > 0) |>
    dplyr::pull(Entrez) |>
    unique()

  # create reference ---------
  create_list = function(binomial_results) {
    list(
      x = binomial_results@statistics@statistic_data$num_positive,
      n = binomial_results@statistics@statistic_data$num_negative +
        binomial_results@statistics@statistic_data$num_positive,
      p = binomial_results@statistics@statistic_data$null_value,
      alternative = binomial_results@statistics@statistic_data$alternative
    )
  }

  # two sided ------------
  binomial_results = binomial_feature_enrichment(
    new(
      "binomial_features",
      positivefc = h_10_pos,
      negativefc = h_10_neg,
      annotation = enrichment_data$bp_annotation
    ),
    p_adjust = "BH"
  )

  input_list = create_list(binomial_results)
  raw_binom = purrr::pmap(input_list, binom.test)
  raw_p = purrr::map_dbl(raw_binom, ~ .x$p.value)
  expect_equal(raw_p, binomial_results@statistics@statistic_data$p)

  expect_snapshot_value(binomial_results, style = "serialize")

  # less than ---------------
  binomial_results = binomial_feature_enrichment(
    new(
      "binomial_features",
      positivefc = h_10_pos,
      negativefc = h_10_neg,
      annotation = enrichment_data$bp_annotation
    ),
    direction = "less",
    p_adjust = "BH"
  )

  input_list = create_list(binomial_results)
  raw_binom = purrr::pmap(input_list, binom.test)
  raw_p = purrr::map_dbl(raw_binom, ~ .x$p.value)
  expect_equal(raw_p, binomial_results@statistics@statistic_data$p)

  expect_snapshot_value(binomial_results, style = "serialize")

  # greater than --------------
  binomial_results = binomial_feature_enrichment(
    new(
      "binomial_features",
      positivefc = h_10_pos,
      negativefc = h_10_neg,
      annotation = enrichment_data$bp_annotation
    ),
    direction = "greater",
    p_adjust = "BH"
  )

  input_list = create_list(binomial_results)
  raw_binom = purrr::pmap(input_list, binom.test)
  raw_p = purrr::map_dbl(raw_binom, ~ .x$p.value)
  expect_equal(raw_p, binomial_results@statistics@statistic_data$p)

  expect_snapshot_value(binomial_results, style = "serialize")
})


test_that("binomial warning works", {
  lipid_features = readRDS("lipid_binomial_testing.rds")
  expect_warning(
    binomial_feature_enrichment(lipid_features, min_features = 6),
    "No annotations had more than"
  )

  expect_no_warning(binomial_feature_enrichment(
    lipid_features,
    min_features = 1
  ))
})

test_that("binomial does normal stuff", {
  h_10_neg = enrichment_data$table10 %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::pull(Entrez) %>%
    unique()

  h_10_pos = enrichment_data$table10 %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::pull(Entrez) %>%
    unique()

  h_48_neg = enrichment_data$table48 |>
    dplyr::filter(logFC < 0) |>
    dplyr::pull(Entrez) |>
    unique()

  h_48_pos = enrichment_data$table48 |>
    dplyr::filter(logFC > 0) |>
    dplyr::pull(Entrez) |>
    unique()
  binomial_h10 = binomial_feature_enrichment(
    new(
      "binomial_features",
      positivefc = h_10_pos,
      negativefc = h_10_neg,
      annotation = enrichment_data$bp_annotation
    ),
    p_adjust = "BH"
  )
  binomial_h48 = binomial_feature_enrichment(
    new(
      "binomial_features",
      positivefc = h_48_pos,
      negativefc = h_48_neg,
      annotation = enrichment_data$bp_annotation
    ),
    p_adjust = "BH"
  )

  # combining ------------
  binomial_comb = combine_enrichments(h10 = binomial_h10, h48 = binomial_h48)
  expect_named(binomial_comb@enriched, c("h10", "h48"))
  n_h10 = sum(grepl("^h10", names(binomial_comb@statistics@statistic_data)))
  n_h48 = sum(grepl("^h48", names(binomial_comb@statistics@statistic_data)))
  expect_equal(n_h10, n_h48)

  # filtering --------------
  binomial_filter = get_significant_annotations(
    binomial_comb,
    padjust <= 0.05
  )

  expect_gt(sum(binomial_filter@statistics@significant@significant[, 1]), 0)
  expect_gt(sum(binomial_filter@statistics@significant@significant[, 2]), 0)
})
