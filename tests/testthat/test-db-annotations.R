test_that("annotations from db work", {
  entrez_db = get_db_annotation()

  expect_equal(
    entrez_db@annotation_features[1],
    list(`GO:0000009` = c("55650", "79087", "85365"))
  )

  ensembl_db = get_db_annotation(feature_type = "ENSEMBL")
  expect_equal(
    ensembl_db@annotation_features[1],
    list(
      `GO:0000009` = c("ENSG00000060642", "ENSG00000182858", "ENSG00000119523")
    )
  )

  entrez_cc = get_db_annotation(annotation_type = "CC")
  cc_desc = head(entrez_cc@description)
  expect_all_true(grepl("complex", cc_desc))

  entrez_pmid = get_db_annotation(annotation_type = "PMID")
  expect_equal(
    entrez_pmid@annotation_features["100103"],
    list(
      `100103` = c(
        "100008587",
        "106632260",
        "106632261",
        "106632262",
        "106632263"
      )
    )
  )
  expect_equal(entrez_pmid@annotation_type, "PMID")
})
