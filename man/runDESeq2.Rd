% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq_functions.R
\name{runDESeq2}
\alias{runDESeq2}
\title{Differential expression analysis with DESeq2}
\usage{
runDESeq2(
  data,
  design = NULL,
  formula = ~1,
  contrasts = NULL,
  prefilter = NULL,
  postfilter = NULL,
  min_counts = 5,
  min_samples = 2,
  ctrlgenes = NULL,
  sizefactors = NULL,
  alpha = 0.05,
  ordered = TRUE,
  df = TRUE,
  ncores = NULL,
  shrink = TRUE,
  ihw = TRUE,
  vst = TRUE,
  rlog = FALSE,
  ...
)
}
\arguments{
\item{data}{Counts matrix}

\item{design}{Experimental design/colData}

\item{formula}{Formula}

\item{contrasts}{Named list of contrasts, specified as c(factor, level, reflevel)}

\item{min_counts}{min_counts in min_samples for filtering}

\item{min_samples}{min_counts in min_samples for filtering}

\item{ctrlgenes}{Control genes for normalization (housekeeping genes)}

\item{sizefactors}{Pre-calculated size factors}

\item{alpha}{Significance level (default = 0.05)}

\item{ordered}{Order results (default = TRUE)}

\item{df}{Return results ass data.frame}

\item{ncores}{Parallel processing}

\item{shrink}{Use ashr for log2FC shrinkage}

\item{ihw}{Use independent hypothesis weighting}

\item{vst}{Add vst-transformed assay}

\item{rlog}{Add rlog-transformed assay}

\item{...}{Parameters passed to DESeq}

\item{filter}{Vector for filtering (subsetting) data}
}
\description{
Differential expression analysis with DESeq2
}
