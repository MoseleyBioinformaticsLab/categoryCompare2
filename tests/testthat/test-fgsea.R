enrichment_data = readRDS("enrichment_data.rds")
bp_annotation = enrichment_data$bp_annotation
h10_ranks = enrichment_data$table10$logFC
names(h10_ranks) = enrichment_data$table10$Entrez

h48_ranks = enrichment_data$table48$logFC
names(h48_ranks) = enrichment_data$table48$Entrez

h48_ranks_sort = sort(h48_ranks)
h48_ranks_order = seq_len(length(h48_ranks_sort))
names(h48_ranks_order) = names(h48_ranks_sort)
mid_point = round(length(h48_ranks_order) / 2)
h48_ranks_order2 = h48_ranks_order - mid_point

test_that("fgsea works", {
  fgsea_features = new(
    "gsea_features",
    ranks = h10_ranks,
    annotation = bp_annotation
  )
  withr::with_seed(1234, {
    fgsea_cc2 = gsea_feature_enrichment(
      fgsea_features,
      min_features = 20,
      max_features = 200
    )
  })
  withr::with_seed(1234, {
    fgsea_fgsea = fgsea::fgsea(
      pathways = bp_annotation@annotation_features,
      stats = h10_ranks,
      minSize = 20,
      maxSize = 200
    )
  })

  withr::with_seed(1234, {
    fgsea_as_fgsea = gsea_feature_enrichment(
      fgsea_features,
      min_features = 20,
      max_features = 200,
      return_type = "fgsea"
    )
  })

  tmp_p = fgsea_cc2@statistics@statistic_data[["p"]]
  names(tmp_p) = NULL
  tmp_paths = fgsea_cc2@statistics@annotation_id
  expect_equal(tmp_paths, fgsea_fgsea[["pathway"]])

  expect_equal(tmp_p, fgsea_fgsea[["pval"]])
  expect_equal(fgsea_fgsea, fgsea_as_fgsea)

  fgsea_from_cc2 = gsea_result_to_fgsea(fgsea_cc2)
  expect_equal(
    fgsea_from_cc2$results |> dplyr::select(-pathway),
    fgsea_as_fgsea |> dplyr::select(-pathway)
  )
  expect_equal(fgsea_from_cc2$ranks, h10_ranks)
})

withr::with_seed(1234, {
  h10_gsea = gsea_feature_enrichment(
    new(
      "gsea_features",
      ranks = h10_ranks,
      annotation = bp_annotation
    ),
    min_features = 20,
    max_features = 200
  )
})

withr::with_seed(1234, {
  h48_gsea = gsea_feature_enrichment(
    new(
      "gsea_features",
      ranks = h48_ranks,
      annotation = bp_annotation
    ),
    min_features = 20,
    max_features = 200
  )
})

withr::with_seed(1234, {
  h48_gsea_order = gsea_feature_enrichment(
    new(
      "gsea_features",
      ranks = h48_ranks_order2,
      annotation = bp_annotation
    ),
    min_features = 20,
    max_features = 200
  )
})

gsea_comb = combine_enrichments(h10 = h10_gsea, h48 = h48_gsea)

test_that("gsea combines", {
  expect_named(gsea_comb@enriched, c("h10", "h48"))
  n_h10 = sum(grepl("^h10", names(gsea_comb@statistics@statistic_data)))
  n_h48 = sum(grepl("^h48", names(gsea_comb@statistics@statistic_data)))
  expect_equal(n_h10, n_h48)
})

test_that("gsea filters", {
  gsea_filter = get_significant_annotations(gsea_comb, padjust <= 0.05)

  expect_gt(sum(gsea_filter@statistics@significant@significant[, 1]), 0)
  expect_gt(sum(gsea_filter@statistics@significant@significant[, 2]), 0)
})
