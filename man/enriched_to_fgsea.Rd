% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gsea.R
\name{enriched_to_fgsea}
\alias{enriched_to_fgsea}
\title{convert enriched object}
\usage{
enriched_to_fgsea(in_enriched)
}
\arguments{
\item{in_enriched}{the enrichment object}
}
\value{
data.table
}
\description{
Takes an `enriched_result`, and converts it to the table expected by `fgsea`. 
This should only be done on those that have `gsea` as the *Enrichment Method*.
}
