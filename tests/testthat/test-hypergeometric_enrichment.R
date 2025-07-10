enrichment_data = readRDS("enrichment_data.rds")

h_10_genes = enrichment_data$table10 %>%
  dplyr::filter(adj.P.Val <= 0.05) %>%
  dplyr::pull(Entrez) |>
  unique()

h_48_genes = enrichment_data$table48 %>%
  dplyr::filter(adj.P.Val <= 0.05) %>%
  dplyr::pull(Entrez) |>
  unique()

universe_genes = enrichment_data$universe

test_that("hypergeometric enrichment stays the same", {
  hypergeom_results = hypergeometric_feature_enrichment(
    new(
      "hypergeometric_features",
      significant = h_10_genes,
      universe = enrichment_data$universe,
      annotation = enrichment_data$bp_annotation
    ),
    p_adjust = "BH"
  )

  expect_snapshot_value(hypergeom_results, style = "serialize")
})

test_that("hypergeometric works", {
  h10_features = new(
    "hypergeometric_features",
    significant = h_10_genes,
    universe = universe_genes,
    annotation = enrichment_data$bp_annotation
  )

  h10_enrich = hypergeometric_feature_enrichment(h10_features)
  expect_equal(
    c("GO:0000002" = 0.6726781),
    h10_enrich@statistics@statistic_data$p[1]
  )
})

h10_hyper = hypergeometric_feature_enrichment(
  new(
    "hypergeometric_features",
    significant = h_10_genes,
    universe = universe_genes,
    annotation = enrichment_data$bp_annotation
  )
)

h48_hyper = hypergeometric_feature_enrichment(
  new(
    "hypergeometric_features",
    significant = h_48_genes,
    universe = universe_genes,
    annotation = enrichment_data$bp_annotation
  )
)

hyper_comb = combine_enrichments(h10 = h10_hyper, h48 = h48_hyper)

test_that("hyper combines", {
  expect_named(hyper_comb@enriched, c("h10", "h48"))
  n_h10 = sum(grepl("^h10", names(hyper_comb@statistics@statistic_data)))
  n_h48 = sum(grepl("^h48", names(hyper_comb@statistics@statistic_data)))
  expect_equal(n_h10, n_h48)
})

test_that("hyper filters", {
  hyper_filter = get_significant_annotations(
    hyper_comb,
    padjust <= 0.05,
    counts > 5
  )

  expect_gt(sum(hyper_filter@statistics@significant@significant[, 1]), 0)
  expect_gt(sum(hyper_filter@statistics@significant@significant[, 2]), 0)
})
