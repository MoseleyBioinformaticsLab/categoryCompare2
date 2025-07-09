run_fun = function() {
  library(categoryCompare2)
  test_stat <- new(
    "statistical_results",
    annotation_id = c("a1", "a2", "a3"),
    statistic_data = list(
      pvalues = c(a1 = 0.01, a2 = 0.5, a3 = 0.0001),
      counts = c(a1 = 5, a2 = 10, a3 = 1),
      odds = c(a1 = 20, a2 = 100, a3 = 0)
    )
  )
  p_cut = 0.05
  eval(p_cut)
  get_significant_annotations(test_stat, pvalues < p_cut)
}

callr::r(run_fun)
