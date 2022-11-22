








#' Import counts from nf-core/rnaseq pipeline
#'
#' @param nfdir Results directory
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
nf_importTX <- function(nfdir, ...){

  cat(crayon::red("Use '$counts' for 3'-tagged sequencing data\n"))

  files <- list.files(path = nfdir, pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
  names(files) <- sapply(strsplit(files, "/", fixed = TRUE), function(x) x[[length(x)-1]])

  tx2gene <- read.delim(file.path(nfdir, "star_salmon", "salmon_tx2gene.tsv"), header = FALSE)
  colnames(tx2gene) = c("tx", "gene_id", "gene_name")

  txi <- tximport::tximport(files, type = "salmon", tx2gene = tx2gene[,c(1,3)], ...)

  txi$counts <- matrix(as.integer(round(txi$counts)), nrow = nrow(txi$counts), dimnames = dimnames(txi$counts))

  txi
}



#' Import samplesheet from nf-core/rnaseq pipeline
#'
#' @param file
#' @param exclude
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
nf_getSamplesheet <- function(file, exclude = c("fastq_1", "fastq_2", "strandedness"), ...){

  df <- read.csv(file, ...)
  df <- df[, !colnames(df) %in% exclude, drop = FALSE]
  tibble::column_to_rownames(df, "sample")

}



#' Get version of the used nf-core/rnaseq pipeline
#'
#' @param ymlfile
#'
#' @return
#' @export
#'
#' @examples
nf_getVersion <- function(ymlfile){
  yml <- yaml::read_yaml(file)
  yml$Workflow$`nf-core/rnaseq`
}




nf_preprocessing <- function(multiqc_dir, design = NULL, ignore = FALSE){

  multiqc_files <- list.files(multiqc_dir, pattern = "multiqc_", full.names = TRUE)
  multiqc_files <- multiqc_files[grep("txt", multiqc_files)]
  names(multiqc_files) <- gsub("multiqc_|\\.txt", "", basename(multiqc_files))
  multiqc_data <- lapply(multiqc_files, function(tmpfile) read.delim(file = tmpfile) )

  multiqc <- data.frame(row.names = multiqc_data$general_stats$Sample, sample = multiqc_data$general_stats$Sample)


  ### Preprocessing steps

  # total reads from fqstqc
  fastqc_raw <- multiqc_data$fastqc[match(multiqc$sample, multiqc_data$fastqc$Sample),]
  multiqc$raw <- fastqc_raw$Total.Sequences

  # trimmed reads from fastqc
  fastqc_trimmed <- multiqc_data$fastqc_1[match(multiqc$sample, multiqc_data$fastqc_1$Sample),]
  multiqc$trimmed <- multiqc$raw - fastqc_trimmed$Total.Sequences

  # mapped and unmapped reads from star
  star <- multiqc_data$star[match(multiqc$sample, multiqc_data$star$Sample),]
  if (!ignore) stopifnot(all(star$total_reads == multiqc$raw - multiqc$trimmed))
  multiqc$unmapped <- rowSums(star[,grepl("unmapped_", colnames(star)) & !grepl("percent", colnames(star))])
  multiqc$multimapped <- star$multimapped + star$multimapped_toomany
  multiqc$mapped <- star$uniquely_mapped
  if (!ignore) stopifnot(all(star$total_reads == multiqc$unmapped + multiqc$multimapped + multiqc$mapped))


  # RSEQC counts tags rather than reads (with tags being mapped fragments of reads, so #tags > #reads)
  # but the total reads also don't add up, so we estimate the fractions
  rseqc <- multiqc_data$rseqc_read_distribution[match(multiqc$sample, multiqc_data$rseqc_read_distribution$Sample),]
  # rseqc$total_reads == star$uniquely_mapped # not matching
  exon_fraction <- (rseqc$cds_exons_tag_count + rseqc$X5_utr_exons_tag_count + rseqc$X3_utr_exons_tag_count)/rseqc$total_tags
  multiqc$nonexon <- star$uniquely_mapped * (1 - exon_fraction)
  multiqc$exon <- star$uniquely_mapped * exon_fraction


  if (!ignore){ with(multiqc,
                     stopifnot(all(nonexon + exon + multimapped + unmapped + trimmed == raw))) }


  ggdf <- multiqc %>% dplyr::select(-c(raw, mapped)) %>% tidyr::pivot_longer(cols = -sample, names_to = "fate", values_to = "reads") %>% as.data.frame()
  ggdf$fate <- factor(ggdf$fate, ordered = TRUE, levels = c("unmapped", "multimapped", "trimmed", "nonexon", "exon"))

  if (!is.null(design)){
    ggdf <- subset(ggdf, sample %in% rownames(design))
    ggdf$sample <- factor(ggdf$sample, ordered = TRUE, levels = rownames(design))
  }

  # plot
  gg <- ggplot2::ggplot(data = ggdf, ggplot2::aes(x = sample , y = reads, fill = fate)) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
          axis.line = ggplot2::element_line(colour = "black"),
          axis.text = ggplot2::element_text(colour = "black"),
          axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::geom_bar(stat="identity", width = 0.65) +
    ggplot2::scale_fill_manual(values = c("#e34c39", "#d99b27", "#e3e649", "#1ea9d4", "#1136ed")) +
    ggplot2::labs(x="", y = "million reads", fill="") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0.02)))

  list(data = multiqc, plot = gg)
}




#' Import and plot preprocessing statistics from nf-core/rnaseq pipeline
#'
#' @param nfdir
#' @param design
#' @param ignore
#'
#' @return
#' @export
#'
#' @examples
nf_summary <- function (nfdir, design = NULL, ignore = FALSE){

  multiqc_dir <- file.path(nfdir, "multiqc", "star_salmon", "multiqc_data")

  multiqc_files <- list.files(multiqc_dir, pattern = "multiqc_",
                              full.names = TRUE)
  multiqc_files <- multiqc_files[grep("txt", multiqc_files)]
  names(multiqc_files) <- gsub("multiqc_|\\.txt", "", basename(multiqc_files))
  multiqc_data <- lapply(multiqc_files, function(tmpfile) read.delim(file = tmpfile))
  multiqc <- data.frame(row.names = multiqc_data$general_stats$Sample,
                        sample = multiqc_data$general_stats$Sample)
  fastqc_raw <- multiqc_data$fastqc[match(multiqc$sample, multiqc_data$fastqc$Sample),
  ]
  multiqc$raw <- fastqc_raw$Total.Sequences

  if (!is.null(multiqc_data$fastqc_1)){
    fastqc_trimmed <- multiqc_data$fastqc_1[match(multiqc$sample,
                                                  multiqc_data$fastqc_1$Sample), ]
    multiqc$trimmed <- multiqc$raw - fastqc_trimmed$Total.Sequences
  }

  star <- multiqc_data$star[match(multiqc$sample, multiqc_data$star$Sample), ]

  if (!ignore)
    stopifnot(all(star$total_reads == multiqc$raw - multiqc$trimmed))
  multiqc$unmapped <- rowSums(star[, grepl("unmapped_", colnames(star)) &
                                     !grepl("percent", colnames(star))])
  multiqc$multimapped <- star$multimapped + star$multimapped_toomany
  multiqc$mapped <- star$uniquely_mapped

  if (!ignore)
    stopifnot(all(star$total_reads == multiqc$unmapped +
                    multiqc$multimapped + multiqc$mapped))
  rseqc <- multiqc_data$rseqc_read_distribution[match(multiqc$sample,
                                                      multiqc_data$rseqc_read_distribution$Sample), ]
  exon_fraction <- (rseqc$cds_exons_tag_count + rseqc$X5_utr_exons_tag_count +
                      rseqc$X3_utr_exons_tag_count)/rseqc$total_tags
  multiqc$nonexon <- star$uniquely_mapped * (1 - exon_fraction)
  multiqc$exon <- star$uniquely_mapped * exon_fraction

  if (!ignore & !is.null(multiqc$trimmed)) {
    with(multiqc, stopifnot(all(nonexon + exon + multimapped + unmapped + trimmed == raw)))
  }

  if (!ignore & is.null(multiqc$trimmed)) {
    with(multiqc, stopifnot(all(nonexon + exon + multimapped + unmapped == raw)))
  }

  ggdf <- multiqc %>% dplyr::select(-c(raw, mapped)) %>% tidyr::pivot_longer(cols = -sample,
                                                                             names_to = "fate", values_to = "reads") %>% as.data.frame()
  ggdf$fate <- factor(ggdf$fate, ordered = TRUE, levels = c("unmapped",
                                                            "multimapped", "trimmed", "nonexon", "exon"))
  if (!is.null(design)) {
    ggdf <- subset(ggdf, sample %in% rownames(design))
    ggdf$sample <- factor(ggdf$sample, ordered = TRUE, levels = rownames(design))
  }
  gg <- ggplot2::ggplot(data = ggdf, ggplot2::aes(x = sample,
                                                  y = reads, fill = fate)) + ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   axis.line = ggplot2::element_line(colour = "black"),
                   axis.text = ggplot2::element_text(colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                       vjust = 0.5)) + ggplot2::geom_bar(stat = "identity",
                                                                                         width = 0.65) + ggplot2::scale_fill_manual(values = c("#e34c39",
                                                                                                                                               "#d99b27", "#e3e649", "#1ea9d4", "#1136ed")) + ggplot2::labs(x = "",
                                                                                                                                                                                                            y = "million reads", fill = "") + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,
                                                                                                                                                                                                                                                                                                               0.02)))
  list(data = multiqc, plot = gg)
}










#' Differential expression analysis with DESeq2
#'
#' @param data
#' @param design
#' @param formula
#' @param contrasts
#' @param filter
#' @param alpha
#' @param ordered
#' @param df
#' @param ncores
#' @param shrink
#' @param ihw
#' @param vst
#' @param rlog
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runDESeq2 <- function(data, design = NULL, formula = ~ 1, contrasts = NULL, filter = NULL,
                      ctrlgenes = NULL, sizefactors = NULL,
                      alpha = 0.05, ordered = TRUE, df = TRUE, ncores = NULL,
                      shrink = TRUE, ihw = TRUE, vst = TRUE, rlog = FALSE, ...){

  stopifnot(requireNamespace("DESeq2", quietly = TRUE))
  stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
  if (ifelse(is.null(ncores), TRUE, ncores > 1)) stopifnot(requireNamespace("BiocParallel", quietly = TRUE))
  if (ihw == TRUE) stopifnot(requireNamespace("IHW", quietly = TRUE))
  if (shrink == TRUE) stopifnot(requireNamespace("ashr", quietly = TRUE))

  results <- list()


  # Parallel setup ----
  if (is.null(ncores)) ncores <- min(c(4, max(c(1, length(contrasts)))))
  bppar <- NULL
  if (ncores > 1){
    BiocParallel::register(BiocParallel::MulticoreParam(workers = ncores))
    bppar <- BiocParallel::bpparam()
    message(paste0("Parallel on ", ncores, " cores..."))
  }


  # Input data ----
  if ("SummarizedExperiment" %in% class(data)){
    if (is.null(design)) design <- data.frame(row.names = colnames(assays(data)[[1]]))
    data <- data[,rownames(design)]
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSet(data, colData = design, design = formula)

  } else if ("matrix" %in% class(data)){
    if (is.null(design)) design <- data.frame(row.names = colnames(data))
    data <- data[,rownames(design)]
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSetFromMatrix(data, colData = design, design = formula)

  } else if ("list" %in% class(data)){
    if (is.null(design)) design <- data.frame(row.names = colnames(data$abundance))
    for (i in seq_along(data)){
      if (!is.null(ncol(data[[i]]))){
        data[[i]] <- data[[i]][,rownames(design)]
      }
    }
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSetFromTximport(data, colData = design, design = formula)

  } else {
    stop("Wrong input data format!")
  }


  # Pre-filtering ----

  if (is.null(filter)){
    dds <- dds[rowSums(DESeq2::counts(dds) > 0, na.rm = TRUE) != 0,]
  } else if (length(filter) == nrow(dds)){
    dds <- dds[naf(filter),]
  } else {
    stop("Index for filtering is of wrong length/format!")
  }


  # Model fitting ----
  if (!is.null(sizefactors)){
    dds$sizeFactor <- sizefactors
  } else if (!is.null(ctrlgenes)){
    if (!all(ctrlgenes %in% rownames(dds))){
      stop("Not all ctrlgenes present in data!")
    }
    dds <- DESeq2::estimateSizeFactors(dds, controlGenes = rownames(dds) %in% ctrlgenes)
  }
  dds <- DESeq2::DESeq(dds, parallel = (ncores > 1), BPPARAM = bppar, ...)


  # Normalized counts ----
  results$normcounts <- DESeq2::counts(dds, normalized = TRUE)
  if (vst == TRUE) results$vst <- SummarizedExperiment::assays(DESeq2::vst(dds, blind = TRUE))[[1]]
  if (rlog == TRUE) results$rlog <- SummarizedExperiment::assays(DESeq2::rlog(dds, blind = TRUE))[[1]]


  # Contrasts ----
  if (!is.null(contrasts)){

    if (ihw == TRUE) filterFun <- IHW::ihw else filterFun <- NULL

    results$results <- lapply(contrasts, function(tmp){

      if (is.null(filterFun)){
        mle <- DESeq2::results(dds, contrast = tmp, alpha = alpha)
      } else {
        mle <- DESeq2::results(dds, filterFun = filterFun, contrast = tmp, alpha = alpha)
      }


      if (shrink == TRUE){
        mmse <- DESeq2::lfcShrink(dds, contrast = tmp, res = mle, type = "ashr", svalue = TRUE)
        mmse <- mmse[rownames(mle),]
        mle$svalue <- mmse$svalue
        mle$log2FCshrink <- mmse$log2FoldChange
      }

      mle$gene <- rownames(mle)
      mle <- as.data.frame(mle)
      mle <- dplyr::relocate(.data = mle, gene)
      mle <- dplyr::rename(.data = mle, log2FC = log2FoldChange, log2FCse = lfcSE)
      if (ordered == TRUE) mle <- dplyr::arrange(mle, padj)
      mle
    })

  }


  # Results ----
  results$dds <- dds
  results$design <- as.data.frame(SummarizedExperiment::colData(dds))
  results

}










#' Differential expression analysis with LIMMA
#'
#' @param data
#' @param design
#' @param formula
#' @param contrasts
#' @param trend
#' @param robust
#' @param p.adj.method
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runLIMMA <- function(data, design, formula = ~ 1, contrasts = NULL, trend = TRUE, robust = FALSE, p.adj.method = "fdr", normalize = FALSE, norm.method = NULL, ...){

  stopifnot(requireNamespace("limma", quietly = TRUE))

  if (normalize == TRUE){
    data <- limma::normalizeBetweenArrays(data, method = norm.method)
  }


  formula <- update(formula,  ~  0 + .)
  mm <- model.matrix(formula, design)

  contrasts_named <- lapply(contrasts, function(contr){  paste(c(paste0(contr[1], contr[2]), paste0(contr[1], contr[3])), collapse = " - ")})
  args <- c(contrasts_named, list(levels = mm))
  contrasts_limma <- do.call(what = limma::makeContrasts, args = args)

  fit <- limma::lmFit(data, mm)
  cfit <- limma::contrasts.fit(fit = fit, contrasts = contrasts_limma)
  efit <- limma::eBayes(fit = cfit, trend = trend, robust = robust, ...)

  coefnames <- colnames(efit$coefficients)
  results <- lapply(setNames(coefnames, coefnames), function(tmpcoef) limma::topTable(efit, coef = tmpcoef, number = Inf, adjust.method = p.adj.method) )
  results <- lapply(results, function(tmpres) data.frame(row.names = rownames(tmpres),
                                                         "id" = rownames(tmpres),
                                                         "mean" = tmpres$AveExpr,
                                                         "stat" = tmpres$t,
                                                         "log2FC" = tmpres$logFC,
                                                         "pvalue" = tmpres$P.Value,
                                                         "padj" = tmpres$adj.P.Val) )


  results
}















#' Collapse duplicate ids
#'
#' @param data matrix/data.frame
#' @param ids vector of row identifiers
#' @param select_by FUN for calculating stat by which duplicates are selected (applied to each row)
#' @param average_by FUN for calculating averages of duplicate rows
#' @param decreasing sort selecting stat increasingly/decreasingly (default = TRUE)
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
collapse <- function(data, ids, select_by = NULL, average_by = NULL, decreasing = TRUE, ...){

  rownames(data) <- NULL

  stopifnot(nrow(data) == length(ids))
  if (any(is.na(ids))){
    message("Warning: 'NA' ids are removed from results.")
  }

  if (sum(nat(duplicated(ids))) == 0){
    message("No ids are duplicated.")
    return(data)
  }

  if (is.null(select_by) & is.null(average_by)){
    message("No method provided. Returning original data.")
    return(data)
  }

  if (!is.null(select_by) & !is.null(average_by)){
    stop("Error: Use only one of 'select_by', 'average_by'!")
  }

  ix <- 1:nrow(data)



  # select_by ------------------------------------

  # Select duplicate with highest (lowest) stat as returned by FUN

  if (!is.null(select_by)){

    stat <- apply(data, 1, FUN = select_by, ...)
    ix_ordered <- order(stat, decreasing = decreasing)
    ix_keep <- sort(ix[ix_ordered][!duplicated(ids[ix_ordered])])

    data <- data[ix_keep,, drop = FALSE]
    rownames(data) <- ids[ix_keep]
    data <- data[!is.na(ids[ix_keep]),, drop = FALSE]
  }



  # average_by ------------------------------------

  # Average duplicates based on FUN

  if (!is.null(average_by)){
    data <- data.frame(data, ids = ids)
    data <- data %>% dplyr::group_by(ids) %>%
      dplyr::summarise(dplyr::across(.fns = average_by, ...))
    data <- as.data.frame(data)
    data <- data[!is.na(data$ids),, drop = FALSE]
    rownames(data) <- data$ids

    data$ids <- NULL
  }


  data[na.omit(match(ids, rownames(data))),, drop = FALSE] # bring into original order


}










#' Get feature ranks from differential abundance testing
#'
#' @param data
#' @param rank_by
#' @param type
#'
#' @return
#' @export
#'
#' @examples
getRanks <- function(data, rank_by = 1, type = NULL){

  # Column 1 is always the main ranking, positive and negative

  rank_by <- rlang::enquo(rank_by)
  rankdf <- dplyr::select(data, !!rank_by)

  if (!is.null(type)){
    stopifnot(length(type) == ncol(rankdf))
    for (i in seq(type)){
      if (type[i] %in% c("+", "positive")) rankdf[[i]] <- rankdf[[i]] * sign(rankdf[[1]])
      if (type[i] %in% c("p", "probability")) rankdf[[i]] <- -log10(rankdf[[i]]) * sign(rankdf[[1]])
    }
  }

  df <- dplyr::arrange_all(rankdf)
  df <- df[rev(1:nrow(df)),,drop = FALSE]

  df$id <- as.character(df[[1]])
  df$id[is.na(df$id)] <- "na"
  df$dup <- df$id %in% df$id[duplicated(df$id)]
  df$newstat <- df[[1]]

  for (id in unique(subset(df, dup == TRUE)$id)){

    tmp <- df[df$id == id,,drop = FALSE]
    ix <- range(which(id == df$id))

    tmp2 <- df[c(ix[1]-1, ix[2]+1),,drop = FALSE]
    maxdiff <- min(abs(tmp2[[1]] - mean(tmp[[1]]))) * 0.4

    tmpranks <- rev(1:nrow(tmp))
    tmpranks <- (tmpranks - median(tmpranks)) / max(abs(tmpranks))

    tmp$newstat <- tmp[[1]] + maxdiff * tmpranks
    df[df$id == id,]$newstat <- tmp$newstat
  }

  ranks <- setNames(df$newstat, rownames(df))
  ranks
}






#' GSEA with clusterProfiler
#'
#' @param data
#' @param genesets
#' @param rank_by
#' @param type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runGSEA <- function(data, genesets = NULL, rank_by = c(stat, baseMean, svalue), type = c("", "+", "p"), as.df = TRUE, ...){

  stopifnot(requireNamespace("clusterProfiler"))
  stopifnot(requireNamespace("fgsea"))

  rank_by <- rlang::enquo(rank_by)

  if (is.null(genesets)){
    genesets <- NA
  }

  if (!"data.frame" %in% class(genesets)){
    genesets <- convertGeneSets(genesets, to = "data.frame")
  }

  ranks <- getRanks(data, rank_by = !!rank_by, type = type)
  ranks <- sort(ranks, decreasing = TRUE)

  results <- clusterProfiler::GSEA(ranks,
                                   seed = 123,
                                   eps = 0,
                                   minGSSize = 1,
                                   maxGSSize = 2000,
                                   TERM2GENE = genesets,
                                   pvalueCutoff = 1,
                                   pAdjustMethod = "fdr",
                                   ...)

  if (as.df == TRUE){
    results <- as.data.frame(results)
    results <- dplyr::select(results, -Description)
    results <- dplyr::rename(results, term = ID, ES = enrichmentScore, padj = p.adjust, qval = qvalues, size = setSize)
  }

  results
}




#' GSEA with fGSEA
#'
#' @param data
#' @param genesets
#' @param rank_by
#' @param type
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
runfGSEA <- function(data, genesets, rank_by = c(stat, baseMean, svalue), type = c("", "+", "p"), seed = 123){

  stopifnot(requireNamespace("fgsea"))

  if ("data.frame" %in% class(genesets)){
    genesets <- convertGeneSets(genesets, to = "list")
  }

  rank_by <- rlang::enquo(rank_by)

  data <- as.data.frame(data)
  ranks <- getRanks(data = data, rank_by = !!rank_by, type = type)

  set.seed(seed)
  results <- fgsea::fgsea(genesets, sort(ranks), eps = 0)
  results <- as.data.frame(results)
  results <- dplyr::rename(results, term = pathway)
  results <- dplyr::arrange(results, pval)
  results
}





#' Single-sample gene set enrichment with GSVA
#'
#' @param data
#' @param genesets
#' @param method
#' @param kcdf
#' @param ncores
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runGSVA <- function(data, genesets = NULL, method = "gsva", kcdf = "Gaussian", ncores = 1, ...){

  stopifnot(requireNamespace("GSVA"))

  if ("data.frame" %in% class(genesets)){
    genesets <- convertGeneSets(genesets, to = "list")
  }

  ncores <- min(8, ncores, parallel::detectCores())
  data <- data.matrix(data)

  results <- GSVA::gsva(data,
                        genesets,
                        method = method,
                        kcdf = kcdf,
                        parallel.sz = ncores,
                        ...)

  results <- data.frame(results)
  results
}





#' Gene set over-representation analysis
#'
#' @param genes vector of genes
#' @param genesets gene set dataframe
#' @param p_cutoff
#' @param ...
#'
#' @return
#' @export
#'
#' @examples getMSigDB("H")$HALLMARK_GLYCOLYSIS |> runORA()
runORA <- function(genes, genesets = NULL, p_cutoff = 0.25, ...){

  stopifnot(requireNamespace("clusterProfiler"))

  if (is.null(genesets)) genesets <- getGOgenes(...)

  results <- clusterProfiler::enricher(gene = genes, TERM2GENE = genesets, pAdjustMethod = "fdr", pvalueCutoff = 1, qvalueCutoff = 1)
  subset(as.data.frame(results), p.adjust <= p_cutoff)
}





getTestData <- function(n = 12, nrow = 100){

  data <- rmat(nrow = nrow, ncol = n, FUN = rpois, lambda = 100)
  design <- data.frame(row.names = colnames(data), group = as.character(partition(1:n, 3)), cov1 = runif(n))

  ix <- 1:round(nrow(data) * 0.1)
  down <- rownames(data)[ix]
  data[ix,design$group == 1] <- data[ix,design$group == 1] * 10

  ix <- round(nrow(data) * 0.9):nrow(data)
  up <- rownames(data)[ix]
  data[ix,design$group == 2] <- data[ix,design$group == 1] * 10


  contrasts <- list(group2vs1 = c("group", 2, 1))

  list(data = data, design = design, contrasts = contrasts, up = up, down = down)


}



getTestGenesets <- function(genes, n = 10){
  genesets <- lapply(1:n, function(i) sample(genes, size = sample(5:30, 1)) )
  names(genesets) <- replicate(paste(c(sample(LETTERS, 1), sample(letters, 5)), collapse = ""), n = n)
  genesets
}









