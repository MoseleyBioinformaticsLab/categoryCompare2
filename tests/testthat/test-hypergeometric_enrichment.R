test_that("hypergeometric enrichment stays the same", {
  enrichment_data = readRDS("enrichment_data.rds")
  
  h_10_genes = enrichment_data$table10 %>%
    dplyr::filter(adj.P.Val <= 0.05) %>%
    dplyr::pull(Entrez)
  
  hypergeom_results = hypergeometric_feature_enrichment(
    new("hypergeom_features", significant = h_10_genes,
        universe = enrichment_data$universe, annotation = enrichment_data$bp_annotation),
    p_adjust = "BH"
  )
  
  expect_snapshot_value(hypergeom_results, style = "serialize")
})


