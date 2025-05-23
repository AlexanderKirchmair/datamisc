% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq_functions.R
\name{get_feature_ranks}
\alias{get_feature_ranks}
\title{Get feature ranking vector from differential abundance testing results}
\usage{
get_feature_ranks(
  data,
  rank_by = "stat",
  rank_type = "",
  id_col = NULL,
  direction_by = NULL,
  ...
)
}
\arguments{
\item{rank_type}{}
}
\description{
Calculate ranking values to sort features from most upregulated to most downregulated for GSEA
}
\examples{
df <- data.frame(id = rgenes(n = 100), stat = rnorm(100), mean = runif(100))
df[1:50,-1] <- df[51:100,-1] # make duplicates
df |> get_feature_ranks(rank_by = c("stat", "mean"), rank_type = c("", ""), id_col="id")

dds <- DESeq2::makeExampleDESeqDataSet(n = 10000, interceptMean = c(2,5))
res <- dds |> runDESeq2(formula = ~ condition, contrasts = list(BvsA = c("condition", "B", "A")))
res$results$BvsA |> get_feature_ranks(rank_by = c("pvalue", "baseMean"), rank_type = c("p", ""), direction_by = "log2FC")
}
