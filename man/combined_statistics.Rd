% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_classes.R, R/combine_enrichments.R
\name{combined_statistics}
\alias{combined_statistics}
\title{combined statistics}
\usage{
combined_statistics(
  statistic_data,
  which_enrichment,
  which_statistic,
  annotation_id,
  significant = NULL,
  measured = NULL,
  use_names = NULL
)

combined_statistics(
  statistic_data,
  which_enrichment,
  which_statistic,
  annotation_id,
  significant = NULL,
  measured = NULL,
  use_names = NULL
)
}
\arguments{
\item{statistic_data}{the data.frame of statistics}

\item{which_enrichment}{which enrichment gave the results}

\item{which_statistic}{which statistics were calculated in each case}

\item{annotation_id}{the annotations for which we are returning statistics}

\item{significant}{the significant annotations}

\item{measured}{the measured annotations}

\item{use_names}{the order of naming}
}
\value{
combined_statistics
}
\description{
holds the results of extracting a bunch of statistics from a \code{\link{combined_enrichment}}
into one entity. This is useful because we want to enable multiple data representations and
simple filtering on the actual \code{data.frame} of statistics, and this provides flexibility
to enable that.

constructor function for the combined_statistics object, makes sure that
empty things get initialized correctly
}
\section{Slots}{

\describe{
\item{\code{statistic_data}}{a \code{data.frame} of all of the statistics from all of the enrichments}

\item{\code{significant}}{a \code{\link{significant_annotations}} object, that may be empty}

\item{\code{which_enrichment}}{a \code{vector} giving which enrichment each column of the statistics came from}

\item{\code{which_statistic}}{a \code{vector} providing which statistic each column contains}
}}

