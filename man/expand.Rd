% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq_functions.R
\name{expand}
\alias{expand}
\title{Expand (demultiplex) rows with multiple ids}
\usage{
expand(data, id_split)
}
\arguments{
\item{data}{Data.frame}

\item{id_split}{List of ID vectors}
}
\description{
Expand (demultiplex) rows with multiple ids
}
\examples{
x <- rmat()
rownames(x) <- paste0(rownames(x), ";s", 1:3)
expand(x, strsplit(rownames(x), split = ";"))
}
