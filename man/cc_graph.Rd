% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_classes.R
\name{cc_graph}
\alias{cc_graph}
\alias{show,cc_graph-method}
\title{cc_graph}
\usage{
cc_graph(graph, significant)

\S4method{show}{cc_graph}(object)

cc_graph(graph, significant)
}
\arguments{
\item{graph}{the \code{\link[graph]{graphNEL}}}

\item{significant}{a matrix indicating which nodes are significant in which experiment}

\item{object}{the cc_graph to show}
}
\description{
A \code{cc_graph} class is a \code{graphNEL} with the added slot of
\code{significant}, a matrix of rows (nodes / annotations) and whether
they were found to be significant in a given enrichment (columns). This
matrix is used for classifying the annotations into different groups, and
generating either pie-charts or coloring the nodes in a visualization.

constructs a \emph{cc_graph} given a \code{\link[graph]{graphNEL}} and a \emph{significant} matrix.
}
\section{Slots}{

\describe{
\item{\code{significant}}{numeric matrix of ones and zeros}
}}

