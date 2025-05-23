% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_functions.R
\name{rgenes}
\alias{rgenes}
\title{Get random genes}
\usage{
rgenes(
  n = 20,
  keytype = "SYMBOL",
  db = org.Hs.eg.db::org.Hs.eg.db,
  type = "protein-coding",
  replace = FALSE,
  prob = NULL,
  ...
)
}
\arguments{
\item{...}{}
}
\description{
Get random genes
}
\examples{
rgenes(100)
}
