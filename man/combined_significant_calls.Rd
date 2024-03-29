% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/significant_results_calls.R
\name{combined_significant_calls}
\alias{combined_significant_calls}
\title{get significant annotations calls}
\usage{
combined_significant_calls(in_results, queries)
}
\arguments{
\item{in_results}{a \code{\link{combined_enrichment}} object}

\item{queries}{a list of queries that can form a call object}
}
\value{
\code{\link{combined_enrichment}} object
}
\description{
In the case where we have a \code{\link{combined_enrichment}} and we want
to get all of the significant annotations from each of them, and put them
together so we can start doing real meta-analysis.
}
\details{
Note that this function returns the original \code{\link{combined_enrichment}} object with a modified
\code{\link{combined_statistics}} slot where the significant annotations have been added in.
}
