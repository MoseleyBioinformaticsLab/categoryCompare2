test_that("gocats annotation importing works", {
  ancestors_file = system.file(
    "extdata",
    "test_data",
    "ancestors.json.gz",
    package = "categoryCompare2"
  )
  namespace_file = system.file(
    "extdata",
    "test_data",
    "namespace.json.gz",
    package = "categoryCompare2"
  )

  ensembl_keys = AnnotationDbi::keys(
    org.Hs.eg.db::org.Hs.eg.db,
    keytype = "ENSEMBL"
  )
  ensembl_uniprot = suppressMessages(AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = ensembl_keys,
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "UNIPROT")
  ))
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

  without_namespace = gocats_to_annotation(
    ancestors_file,
    namespace_file = NULL
  )
  has_bp = all(grepl("^BP|^CC|^MF", without_namespace@description))
  expect_true(!has_bp)
  with_namespace = gocats_to_annotation(ancestors_file, namespace_file)
  has_namespace = all(grepl("^BP|^CC|^MF", with_namespace@description))
  expect_true(has_namespace)

  expect_equal(without_namespace@annotation_features[[1]][1], "O60313")
  expect_equal(without_namespace@annotation_type, "gocatsGO")
  expect_equal(without_namespace@feature_type, "Uniprot")
  expect_equal(length(without_namespace@annotation_features), 12221)

  without_limits = gocats_to_annotation(
    ancestors_file,
    namespace_file = NULL,
    feature_min = 1,
    feature_max = Inf
  )
  expect_equal(length(without_limits@annotation_features), 22415)

  with_translation = gocats_to_annotation(
    ancestors_file,
    namespace_file,
    feature_type = "ENSEMBL",
    annotation_type = "whatever",
    feature_translation = ensembl_uniprot
  )
  expect_equal(with_translation@annotation_features[[1]][1], "ENSG00000143799")
  expect_equal(with_translation@annotation_type, "whatever")
  expect_equal(with_translation@feature_type, "ENSEMBL")
})
