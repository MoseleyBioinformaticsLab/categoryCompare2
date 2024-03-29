% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_visualization.R
\name{assign_colors}
\alias{assign_colors}
\title{assign colors}
\usage{
assign_colors(in_assign, type = "experiment")
}
\arguments{
\item{in_assign}{the \code{\link{node_assign}} object generated from a \code{\link{cc_graph}}}

\item{type}{either "group" or "experiment"}
}
\value{
node_assign with colors
}
\description{
given a \code{\link{node_assign}}, assign colors to either the independent groups
of unique annotations, or to each of the experiments independently.
}
