% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq_functions.R
\name{runTEST}
\alias{runTEST}
\title{Differential expression analysis with custom test}
\usage{
runTEST(
  data,
  design,
  contrasts = NULL,
  FUN = t.test,
  pull = c(p = "p.value", stat = "statistic"),
  ...
)
}
\arguments{
\item{...}{}
}
\description{
Differential expression analysis with custom test
}
