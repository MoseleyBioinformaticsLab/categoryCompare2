% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotations.R
\name{annotation_gene_table}
\alias{annotation_gene_table}
\title{annotation to genes}
\usage{
annotation_gene_table(
  combined_enrichment,
  annotations = NULL,
  use_db = NULL,
  input_type = "ENTREZID",
  gene_info = c("SYMBOL", "GENENAME")
)
}
\arguments{
\item{combined_enrichment}{combined enrichment object}

\item{annotations}{which annotations to grab features from}

\item{use_db}{the annotation database}

\item{input_type}{what type of gene id was it?}

\item{gene_info}{what type of info to return for each gene}
}
\value{
data.frame
}
\description{
Creates a tabular output of annotations to genes providing lookup of which
genes are contributing to a particular annotation.
}
