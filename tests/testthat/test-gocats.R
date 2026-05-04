test_that("gocats annotation importing works", {
  ensembl_uniprot = enrichment_data$uniprot_to_ensembl
  names(ensembl_uniprot) = c("to", "from")
  ensembl_uniprot = ensembl_uniprot[
    !(is.na(ensembl_uniprot$to)) & !is.na(ensembl_uniprot$from),
  ]
  bad_translation = ensembl_uniprot
  names(bad_translation) = c("other", "one")
  list_translation = as.list(ensembl_uniprot)

  expect_error(gocats_to_annotation("random_file"))
  expect_warning(
    gocats_to_annotation(ancestors_file, "random_file"),
    "random_file does not exist"
  )
  expect_error(
    gocats_to_annotation(
      ancestors_file,
      feature_translation = list_translation
    ),
    "must be a data.frame"
  )
  expect_error(
    gocats_to_annotation(ancestors_file, feature_translation = bad_translation),
    "must contain the columns"
  )

  without_namespace_description = suppressMessages(gocats_to_annotation(
    ancestors_file,
    namespace_file = NULL,
    add_description = "no"
  ))

  expect_equal(
    without_namespace_description@annotation_features[[1]][1],
    "O60313"
  )
  expect_equal(without_namespace_description@annotation_type, "gocatsGO")
  expect_equal(without_namespace_description@feature_type, "Uniprot")
  # 11927 updated for Bioc 3.22 comes from latest version of org.Hs.eg.db, will fail locally
  expect_equal(length(without_namespace_description@annotation_features), 12221)

  without_namespace = gocats_to_annotation(
    ancestors_file,
    namespace_file = NULL
  )
  has_nonamespace = grepl("^BP\\:|^CC\\:|^MF\\:", without_namespace@description)
  expect_true(all(!has_nonamespace))

  with_namespace = gocats_to_annotation(ancestors_file, namespace_file)
  has_namespace = all(grepl("^BP\\:|^CC\\:|^MF\\:", with_namespace@description))
  expect_true(has_namespace)

  without_limits = gocats_to_annotation(
    ancestors_file,
    namespace_file = NULL,
    feature_min = 1,
    feature_max = Inf
  )
  # 21633 for Bioc 3.22 comes from latest version of org.Hs.eg.db, this fails locally
  expect_gt(
    length(without_limits@description),
    length(without_namespace@description)
  )

  with_translation = suppressMessages(gocats_to_annotation(
    ancestors_file,
    namespace_file,
    feature_type = "ENSEMBL",
    annotation_type = "whatever",
    feature_translation = ensembl_uniprot,
    add_description = "no"
  ))
  expect_equal(with_translation@annotation_features[[1]][1], "ENSG00000143799")
  expect_equal(with_translation@annotation_type, "whatever")
  expect_equal(with_translation@feature_type, "ENSEMBL")
})

test_that("json exporting and importing works properly", {
  without_namespace = gocats_to_annotation(
    ancestors_file,
    namespace_file = NULL
  )
  withr::with_file("test_all.json", {
    annotation_2_json(without_namespace, json_file = "test_all.json")
    in_annotation = json_2_annotation("test_all.json")
    expect_equal(without_namespace, in_annotation)
  })

  withr::with_file("test_small.json", {
    without_namespace2 = without_namespace
    without_namespace2@annotation_features[[
      1
    ]] = without_namespace@annotation_features[[1]][1]
    annotation_2_json(without_namespace2, json_file = "test_small.json")
    line_data = base::readLines("test_small.json", 20)
    expect_true(grepl("\\[.*\\]", line_data[3]))

    without_namespace2@counts[1] = 1
    in_annotation = json_2_annotation("test_small.json")
    expect_equal(without_namespace2, in_annotation)
  })
})

test_that("json reversal works", {
  withr::with_file("reversed.json", {
    out_res = json_annotation_reversal(
      ancestors_file,
      out_file = "reversed.json"
    )
    tmp_json = jsonlite::fromJSON("reversed.json", simplifyVector = FALSE)
    expect_equal(
      names(tmp_json),
      c("annotation_features", "counts", "annotation_type", "feature_type")
    )
    expect_equal(
      head(names(tmp_json$annotation_features)),
      c(
        "GO:0000002",
        "GO:0000003",
        "GO:0000009",
        "GO:0000010",
        "GO:0000012",
        "GO:0000014"
      )
    )
  })
})
