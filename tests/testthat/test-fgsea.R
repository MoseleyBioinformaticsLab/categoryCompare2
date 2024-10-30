test_that("fgsea works", {
  enrichment_data = readRDS("enrichment_data.rds")
  bp_annotation = enrichment_data$bp_annotation
  ranks = enrichment_data$table10$logFC
  names(ranks) = enrichment_data$table10$Entrez
  
  fgsea_features = new("gsea_features",
                        ranks = ranks,
                      annotation = bp_annotation)
  withr::with_seed(1234, {
    fgsea_cc2 = gsea_feature_enrichment(fgsea_features,
      min_features = 20,
    max_features = 200)
  })
  withr::with_seed(1234, {
    fgsea_fgsea = fgsea::fgsea(pathways = bp_annotation@annotation_features,
      stats = ranks,
      minSize = 20,
      maxSize = 200)
  })

  withr::with_seed(1234, {
    fgsea_as_fgsea = gsea_feature_enrichment(fgsea_features,
      min_features = 20,
    max_features = 200,
    return_type = "fgsea")
  })
  
  tmp_p = fgsea_cc2@statistics@statistic_data[["p"]]
  names(tmp_p) = NULL
  tmp_paths = fgsea_cc2@statistics@annotation_id
  expect_equal(tmp_paths, fgsea_fgsea[["pathway"]])
  
  expect_equal(tmp_p,
              fgsea_fgsea[["pval"]])
  expect_equal(fgsea_fgsea, fgsea_as_fgsea)

  fgsea_from_cc2 = enriched_to_fgsea(fgsea_cc2)
  expect_equal(fgsea_from_cc2$results |> dplyr::select(-pathway),
               fgsea_as_fgsea |> dplyr::select(-pathway))
  expect_equal(fgsea_from_cc2$ranks, ranks)
})
