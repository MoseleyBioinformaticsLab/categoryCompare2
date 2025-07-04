% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_enrichments.R
\name{extract_statistics,combined_enrichment-method}
\alias{extract_statistics,combined_enrichment-method}
\title{extract statistics}
\usage{
\S4method{extract_statistics}{combined_enrichment}(in_results)
}
\arguments{
\item{in_results}{the \code{\link{combined_enrichment}} object}
}
\value{
combined_statistics
}
\description{
extract all statistics from a \code{\link{combined_enrichment}} object and
create a \code{\link{combined_statistics}} where each statistic from the underlying
\code{\link{statistical_results}} object in each of the enrichments
is named according to which enrichment it was in and what statistic it was.
}
