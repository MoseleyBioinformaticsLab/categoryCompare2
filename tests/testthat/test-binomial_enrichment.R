enrichment_data = readRDS("enrichment_data.rds")

h_10_neg = enrichment_data$table10 %>%
  dplyr::filter(logFC < 0) %>%
  dplyr::pull(Entrez) %>%
  unique()

h_10_pos = enrichment_data$table10 %>%
  dplyr::filter(logFC > 0) %>%
  dplyr::pull(Entrez) %>%
  unique()

create_list = function(binomial_results){
  list(x = binomial_results@statistics@statistic_data$statistic,
       n = binomial_results@statistics@statistic_data$parameter,
       p = binomial_results@statistics@statistic_data$null_value,
       alternative = binomial_results@statistics@statistic_data$alternative)
}

test_that("binomial two sided works", {
  
  binomial_results = binomial_feature_enrichment(
    new("binomial_features", positivefc = h_10_pos,
        negativefc = h_10_neg, 
        annotation = enrichment_data$bp_annotation),
    p_adjust = "BH"
  )
  
  input_list = create_list(binomial_results)
  raw_binom = purrr::pmap(input_list, binom.test)
  raw_p = purrr::map_dbl(raw_binom, ~ .x$p.value)
  expect_equal(raw_p, binomial_results@statistics@statistic_data$p)
  
  expect_snapshot_value(binomial_results, style = "serialize")
})

test_that("binomial less works", {
  
  binomial_results = binomial_feature_enrichment(
    new("binomial_features", positivefc = h_10_pos,
        negativefc = h_10_neg, 
        annotation = enrichment_data$bp_annotation),
    direction = "less",
    p_adjust = "BH"
  )
  
  input_list = create_list(binomial_results)
  raw_binom = purrr::pmap(input_list, binom.test)
  raw_p = purrr::map_dbl(raw_binom, ~ .x$p.value)
  expect_equal(raw_p, binomial_results@statistics@statistic_data$p)
  
  expect_snapshot_value(binomial_results, style = "serialize")
})

test_that("binomial greater works", {
  
  binomial_results = binomial_feature_enrichment(
    new("binomial_features", positivefc = h_10_pos,
        negativefc = h_10_neg, 
        annotation = enrichment_data$bp_annotation),
    direction = "greater",
    p_adjust = "BH"
  )
  
  input_list = create_list(binomial_results)
  raw_binom = purrr::pmap(input_list, binom.test)
  raw_p = purrr::map_dbl(raw_binom, ~ .x$p.value)
  expect_equal(raw_p, binomial_results@statistics@statistic_data$p)
  
  expect_snapshot_value(binomial_results, style = "serialize")
})
