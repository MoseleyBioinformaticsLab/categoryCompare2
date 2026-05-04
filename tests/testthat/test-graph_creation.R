# tests creating the cc_graph objects

#library(graph)

# overlap calculations

test_that("various graph properties are correct", {
  group1 <- letters[1:6]
  group2 <- group1
  group3 <- letters[1:3]
  group4 <- letters[c(1, 2, 8, 10)]

  # jaccard calculations are correct --------
  expect_equal(jaccard_coefficient(group1, group2), 1)
  expect_equal(jaccard_coefficient(group1, group3), 1)
  expect_equal(jaccard_coefficient(group1, group4), 0.5)
  expect_equal(jaccard_coefficient(group2, group3), 1)
  expect_equal(jaccard_coefficient(group2, group4), 0.5)
  expect_equal(jaccard_coefficient(group3, group4), 2 / 3)
  expect_equal(jaccard_coefficient(group1, character()), NaN)

  # overlap calculations are correct --------
  expect_equal(overlap_coefficient(group1, group2), 1)
  expect_equal(overlap_coefficient(group1, group3), 0.5)
  expect_equal(overlap_coefficient(group1, group4), 0.25)
  expect_equal(overlap_coefficient(group2, group3), 0.5)
  expect_equal(overlap_coefficient(group2, group4), 0.25)
  expect_equal(overlap_coefficient(group3, group4), 0.4)
  expect_equal(overlap_coefficient(group1, character()), 0)

  # combined calculations are correct ----------
  expect_equal(combined_coefficient(group1, group2), 1)
  expect_equal(combined_coefficient(group1, group3), 0.75)
  expect_equal(combined_coefficient(group1, group4), 0.375)
  expect_equal(combined_coefficient(group2, group3), 0.75)
  expect_equal(combined_coefficient(group2, group4), 0.375)
  expect_equal(combined_coefficient(group3, group4), 8 / 15)
  expect_equal(combined_coefficient(group1, character()), NaN)

  # generating simple graphs ---------
  group5 <- letters[10:12]
  group6 <- letters[11:12]

  feature_list <- list(
    g1 = group1,
    g2 = group2,
    g3 = group3,
    g4 = group4,
    g5 = group5,
    g6 = group6
  )

  # jaccard graph equality -----------
  jaccard_graph <- generate_annotation_similarity_graph(feature_list, "jaccard")
  expect_equal(graph::numNodes(jaccard_graph), 6)
  expect_equal(graph::numEdges(jaccard_graph), 8)

  test_jaccard = generate_annotation_similarity_graph(feature_list, "jaccard")
  expect_snapshot(test_jaccard@edgeData)

  # overlap graph equality -----------
  test_overlap = generate_annotation_similarity_graph(feature_list, "overlap")
  expect_snapshot(test_overlap@edgeData)

  # combined graph equality ------------
  test_combined = generate_annotation_similarity_graph(feature_list, "combined")
  expect_snapshot(test_combined@edgeData)
})


test_that('extraction of significant to graph works', {
  set.seed(1234)
  random_genes = sample(enrichment_data$universe, 100)
  non_enriched = hypergeometric_feature_enrichment(
    new(
      "hypergeometric_features",
      significant = random_genes,
      universe = enrichment_data$universe,
      annotation = enrichment_data$bp_annotation
    ),
    p_adjust = "BH"
  )

  # checking for warnings -------------
  non_list = combine_enrichments(c1 = non_enriched)
  non_sig = get_significant_annotations(non_list, padjust <= 0.01)
  expect_warning(
    generate_annotation_graph(non_sig, low_cut = 2, hi_cut = 1000),
    regexp = "Nothing significant"
  )

  # and that we get a 27 member network -------------
  non_sig2 = get_significant_annotations(non_list, p <= 0.01)
  sig2_graph = generate_annotation_graph(non_sig2, low_cut = 2, hi_cut = 1000)
  expect_equal(graph::numNodes(sig2_graph), 27)

  # checking that tooltips get added ------------
  sig2_graph_tooltips = add_tooltip(
    sig2_graph,
    description = non_sig2@annotation@description[graph::nodes(sig2_graph)]
  )
  expect_equal(
    graph::nodeData(sig2_graph_tooltips, "GO:0000086", "tooltip"),
    list(
      `GO:0000086` = "GO:0000086\nG2/M transition of mitotic cell cycle\nG2/M transition of mitotic cell cycle"
    )
  )

  # can generate communities ---------------
  sig2_comms = assign_communities(sig2_graph_tooltips)
  sig2_labels = label_communities(sig2_comms, non_sig2@annotation)
  expect_equal(length(sig2_comms), 7)
  expect_equal(length(sig2_labels), 5)
  expect_equal(
    sig2_labels[[1]],
    list(
      rep = "GO:1905114",
      label = "cell surface receptor signaling pathway involved in cell-cell signaling",
      members = c(
        "GO:0002524",
        "GO:0021545",
        "GO:0038166",
        "GO:0042249",
        "GO:1905114"
      )
    )
  )

  bare_table = table_from_graph(sig2_graph_tooltips)
  expect_equal(nrow(bare_table), 27)

  tmp_combinations = annotation_combinations(sig2_graph_tooltips)
  tmp_combinations = assign_colors(tmp_combinations)
  comm_table = table_from_graph(
    sig2_graph_tooltips,
    in_assign = tmp_combinations
  )
  expect_equal(nrow(comm_table), nrow(bare_table))

  group_table = table_from_graph(
    sig2_graph_tooltips,
    tmp_combinations,
    sig2_labels
  )
  expect_gt(nrow(group_table), nrow(bare_table))
  expect_equal(group_table$c1.p[1], NA_real_)
})

test_that("check packages installed works", {
  expect_null(check_package_installed("stats"))
  # withr::with_temp_libpaths(
  #   {
  #     expect_error(check_package_installed("Cairo"), "is not installed")
  #   },
  #   action = "replace"
  # )
})

test_that("colors come out correctly", {
  expect_equal(
    generate_colors(4),
    c("#FF7A9E", "#AAB500", "#00CEB7", "#AC9AFF")
  )
  expect_equal(
    generate_colors(5),
    c("#FF7A9E", "#00CA66", "#DA88FF", "#CAAA00", "#00C6F2")
  )
})

test_that("modifying graphs works", {
  group1 <- letters[1:6]
  group2 <- group1
  group3 <- letters[1:3]
  group4 <- letters[c(1, 2, 8, 10)]

  group5 <- letters[10:12]
  group6 <- letters[11:12]

  feature_list <- list(
    g1 = group1,
    g2 = group2,
    g3 = group3,
    g4 = group4,
    g5 = group5,
    g6 = group6
  )

  test_combined = generate_annotation_similarity_graph(feature_list)
  # removing edges is correct ------------
  # remove edges < 0.75, results in 4 edges left
  remove_07 = remove_edges(test_combined, 0.75)
  expect_equal(length(graph::edgeData(remove_07, attr = "weight")), 4)
})
