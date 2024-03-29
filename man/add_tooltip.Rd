% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_visualization.R
\name{add_tooltip}
\alias{add_tooltip}
\title{add tooltip}
\usage{
add_tooltip(
  in_graph,
  node_data = c("name", "description"),
  description,
  separator = "\\n"
)
}
\arguments{
\item{in_graph}{the graph to work with}

\item{node_data}{which pieces of node data to use}

\item{description}{other descriptive text to use}

\item{separator}{what separator to use for the tooltip}
}
\value{
the graph with a new nodeData member "tooltip"
}
\description{
before passing to Cytoscape, add a tooltip attribute to the graph
}
