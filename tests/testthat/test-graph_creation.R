# tests creating the cc_graph objects

#library(graph)

# overlap calculations
group1 <- letters[1:6]
group2 <- group1
group3 <- letters[1:3]
group4 <- letters[c(1, 2, 8, 10)]


test_that("jaccard calculations are correct", {
  expect_equal(jaccard_coefficient(group1, group2), 1)
  expect_equal(jaccard_coefficient(group1, group3), 1)
  expect_equal(jaccard_coefficient(group1, group4), 0.5)
  expect_equal(jaccard_coefficient(group2, group3), 1)
  expect_equal(jaccard_coefficient(group2, group4), 0.5)
  expect_equal(jaccard_coefficient(group3, group4), 2/3)
  expect_equal(jaccard_coefficient(group1, character()), NaN)
})

test_that("overlap calculations are correct", {
  expect_equal(overlap_coefficient(group1, group2), 1)
  expect_equal(overlap_coefficient(group1, group3), 0.5)
  expect_equal(overlap_coefficient(group1, group4), 0.25)
  expect_equal(overlap_coefficient(group2, group3), 0.5)
  expect_equal(overlap_coefficient(group2, group4), 0.25)
  expect_equal(overlap_coefficient(group3, group4), 0.4)
  expect_equal(overlap_coefficient(group1, character()), 0)
})

test_that("combined calculations are correct", {
  expect_equal(combined_coefficient(group1, group2), 1)
  expect_equal(combined_coefficient(group1, group3), 0.75)
  expect_equal(combined_coefficient(group1, group4), 0.375)
  expect_equal(combined_coefficient(group2, group3), 0.75)
  expect_equal(combined_coefficient(group2, group4), 0.375)
  expect_equal(combined_coefficient(group3, group4), 8/15)
  expect_equal(combined_coefficient(group1, character()), NaN)
})

# generating simple graphs
group5 <- letters[10:12]
group6 <- letters[11:12]

feature_list <- list(g1 = group1, g2 = group2, g3 = group3,
                     g4 = group4, g5 = group5, g6 = group6)

jaccard_graph <- generate_annotation_similarity_graph(feature_list, "jaccard")
test_that("number nodes and edges are correct", {
  expect_equal(graph::numNodes(jaccard_graph), 6)
  expect_equal(graph::numEdges(jaccard_graph), 8)
})

test_that("jaccard graph is correct", {
  test_jaccard = generate_annotation_similarity_graph(feature_list, "jaccard")
  expect_snapshot(test_jaccard@edgeData)
})

test_that("overlap graph is correct", {
  test_overlap = generate_annotation_similarity_graph(feature_list, "overlap")
  expect_snapshot(test_overlap@edgeData)
})

test_that("combined graph is correct", {
  test_combined = generate_annotation_similarity_graph(feature_list, "combined")
  expect_snapshot(test_combined@edgeData)
})
