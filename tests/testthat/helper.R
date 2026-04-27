json_bp_annotation = categoryCompare2::json_2_annotation(system.file(
  "extdata",
  "test_data",
  "bp_annotations.gz",
  package = "categoryCompare2"
))

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

enrichment_data = readRDS("enrichment_data.rds")
