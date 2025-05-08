



#' Import counts from nf-core/rnaseq pipeline using tximport
#'
#' @param nfdir results directory
#' @param samplesheet samplesheet data frame
#' @param subdir star_salmon
#' @param pattern quant.sf (salmon output)
#' @param tx2gene tx2gene.tsv (salmon output)
#' @param type
#' @param gene
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
nf_tximport <- function(nfdir, samplesheet = NULL, subdir = "star_salmon", pattern = "quant.sf", tx2gene = "tx2gene.tsv", type = "salmon", gene = "SYMBOL", ...){

  stopifnot(requireNamespace("tximport", quietly = TRUE))

  cat(crayon::blue("'$counts' is the estimate of the expected number of reads\n"))
  cat(crayon::blue("'$abundance' is the estimate of the relative abundance in units of Transcripts Per Million (TPM)\n"))
  cat(crayon::blue("'$length' is the estimate of the effective transcript length in each sample\n"))
  cat(crayon::red("For 3'-tagged sequencing data, run salmon with '--noLengthCorrection'\n"))

  files <- list.files(path = file.path(nfdir, subdir), pattern = pattern, recursive = TRUE, full.names = TRUE)
  names(files) <- sapply(strsplit(files, "/", fixed = TRUE), function(x) x[[length(x) - 1]])

  tx2gene <- read.delim(file.path(nfdir, subdir, tx2gene), header = FALSE, col.names = c("ENSEMBLTRANS", "ENSEMBL", "SYMBOL"))

  txi <- tximport::tximport(files, type = type, tx2gene = tx2gene[, c("ENSEMBLTRANS", gene)], ...)
  txi$counts <- matrix(as.integer(round(txi$counts)), nrow = nrow(txi$counts), dimnames = dimnames(txi$counts))

  if (!is.null(samplesheet)){
    ix <- !tolower(names(txi)) %in% "countsfromabundance"
    ids <- rownames(samplesheet)
    if (!all( sapply(txi[ix], function(X) setequal(ids, colnames(X))) )){
      stop("Error: Mismatching sample names!")
    }
    txi[ix] <- lapply(txi[ix], function(X) X[,rownames(samplesheet), drop = FALSE] )
  }

  txi
}



#' Get versions of tools used in the nf-core/rnaseq pipeline
#'
#' @param ymlfile
#'
#' @return
#' @export
#'
#' @examples
nf_versions <- function(ymlfile){
  yml <- yaml::read_yaml(file)
  yml
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
#' @param data Counts matrix
#' @param design Experimental design/colData
#' @param formula Formula
#' @param contrasts Named list of contrasts, specified as c(factor, level, reflevel)
#' @param lrt_reduced Reduced formula for LRT test
#' @param prefilter Filter before normalization
#' @param postfilter Filter after normalization
#' @param min_counts min_counts in min_samples for filtering
#' @param min_samples min_counts in min_samples for filtering
#' @param ctrlgenes Control genes for normalization (housekeeping genes)
#' @param sizefactors Pre-calculated size factors
#' @param RUV RUV batch effect correction: 'list(empirical = *genes*, group = *column*, n = *n_vars*)'
#' @param SVA SVA batch effect correction: 'list(reduced = *formula*, n = *n_vars*)'
#' @param alpha Significance level (default = 0.05)
#' @param ordered Order results (default = TRUE)
#' @param df Return results as data.frame
#' @param ncores Parallel processing
#' @param shrink Use ashr for log2FC shrinkage
#' @param ihw Use independent hypothesis weighting
#' @param vst Add vst-transformed assay
#' @param rlog Add rlog-transformed assay
#' @param minReplicatesForReplace Required number of replicates to replace outliers
#' @param fitType Fit type of dispresion estimate (parametric, local, mean or glmGamPoi)
#' @param ... Parameters passed to DESeq
#'
#' @return
#' @export
#'
#' @examples
#' dds <- DESeq2::makeExampleDESeqDataSet(n = 10000, interceptMean = c(2,5))
#' dds |> runDESeq2(formula = ~ condition, contrasts = list(BvsA = c("condition", "B", "A")))
#' dds |> runDESeq2(formula = ~ condition, lrt_reduced = ~ 1)
#' dds |> runDESeq2(formula = ~ condition, contrasts = list(BvsA = c("condition", "B", "A")), vst = TRUE, ncores = 5)
#' dds |> runDESeq2(formula = ~ condition, RUV = list(empirical = sample(rownames(dds), 10)), contrasts = list(BvsA = c("condition", "B", "A")))
#' dds |> runDESeq2(formula = ~ condition, SVA = list(reduced = ~ 1), contrasts = list(BvsA = c("condition", "B", "A")))
runDESeq2 <- function(data, design = NULL, formula = ~ 1, contrasts = NULL, lrt_reduced = NULL,
                      prefilter = NULL, postfilter = NULL, min_counts = 5, min_samples = 2,
                      ctrlgenes = NULL, sizefactors = NULL,
                      RUV = list(), SVA = list(),
                      alpha = 0.05, ordered = TRUE, df = TRUE, ncores = NULL,
                      shrink = TRUE, ihw = TRUE, vst = FALSE, rlog = FALSE,
                      minReplicatesForReplace = 7, fitType = "parametric", ...){

  datamisc::colorcat("DESeq2 differential expression analysis", col = "blue")
  datamisc::colorcat("-use 'prefilter' to filter genes before normalization", col = "blue")
  datamisc::colorcat("-use 'ctrlgenes' for normalization to control genes", col = "blue")
  datamisc::colorcat("-use 'RUV' or 'SVA' for batch effect correction", col = "blue")
  datamisc::colorcat("-use 'sizefactors' to directly pass normalization factors", col = "blue")
  datamisc::colorcat("-use 'postfilter' to filter genes after normalization", col = "blue")
  datamisc::colorcat("-use 'plotDispEsts(dds)' to check dispersion estimates", col = "blue")


  stopifnot(requireNamespace("DESeq2", quietly = TRUE))
  stopifnot(requireNamespace("SummarizedExperiment", quietly = TRUE))
  if (ifelse(is.null(ncores), TRUE, ncores > 1)) stopifnot(requireNamespace("BiocParallel", quietly = TRUE))
  if (ihw == TRUE) stopifnot(requireNamespace("IHW", quietly = TRUE))
  if (shrink == TRUE) stopifnot(requireNamespace("ashr", quietly = TRUE))
  if (length(RUV) > 0) stopifnot(requireNamespace("RUVSeq", quietly = TRUE))
  if (length(SVA) > 0) stopifnot(requireNamespace("sva", quietly = TRUE))

  library(DESeq2)

  results <- list()


  # Parallel setup ----
  if (is.null(ncores)) ncores <- min(c(10, max(c(1, length(contrasts)))))
  bppar <- NULL
  if (ncores > 1){
    if (tolower(.Platform$OS.type) == "windows"){
      BiocParallel::register(BiocParallel::SnowParam(workers = ncores))
    } else {
      BiocParallel::register(BiocParallel::MulticoreParam(workers = ncores))
    }
    bppar <- BiocParallel::bpparam()
    message(paste0("Parallel on ", ncores, " cores..."))

    lapply_fun <- function(X, FUN, BPPARAM = NULL, ...){
      BiocParallel::bplapply(X, FUN = FUN, BPPARAM = BPPARAM, ...)
    }

  } else {
    lapply_fun <- function(X, FUN, BPPARAM = NULL, ...){
      lapply(X, FUN = FUN)
    }
  }


  # Input data ----
  if ("SummarizedExperiment" %in% class(data)){
    cat(crayon::blue("Using 'SummarizedExperiment' as input for 'DESeqDataSet'\n"))
    if (is.null(design)) design <- data.frame(row.names = colnames(assays(data)[[1]]))
    data <- data[,rownames(design)]
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSet(data, colData = design, design = formula)

  } else if (any(class(data) %in% c("matrix", "data.frame"))){
    cat(crayon::blue("Using raw counts matrix as input for 'DESeqDataSetFromMatrix'\n"))
    if (is.null(design)) design <- data.frame(row.names = colnames(data))
    data <- data[,rownames(design)]
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSetFromMatrix(data, colData = design, design = formula)

  } else if ("list" %in% class(data)){
    cat(crayon::blue("Using tximport list as input for 'DESeqDataSetFromTximport'\n"))
    if (is.null(design)) design <- data.frame(row.names = colnames(data$abundance))
    for (i in seq_along(data)){
      if (!is.null(ncol(data[[i]]))){
        data[[i]] <- data[[i]][,rownames(design)]
      }
    }
    design <- droplevels(design)
    dds <- DESeq2::DESeqDataSetFromTximport(data, colData = design, design = formula)

  } else if ("DESeqDataSet" %in% class(data)){
    cat(crayon::blue("Using a 'DESeqDataSet' object as input\n"))
      design <- SummarizedExperiment::colData(data)
      dds <- data

  } else {
    stop("Wrong input data format!")
  }


  # Pre-filtering ----
  if (!is.null(prefilter)){
    if (length(prefilter) == nrow(dds) & is.logical(prefilter)){
      dds <- dds[naf(prefilter),]
    } else if (is.character(prefilter)){
      dds <- dds[intersect(prefilter, rownames(dds)),]
    } else {
      stop("Index for pre-filtering is of wrong length/format!")
      }
  }

  if (!is.null(min_counts) & !is.null(min_samples)){
    dds <- dds[rowSums(DESeq2::counts(dds) > min_counts, na.rm = TRUE) >= min_samples,]
  }


  # RUVSeq ----
  if (length(RUV) > 0){
    if (!class(RUV) == "list") stop("Error: 'RUV' must be a named list.")
    if ("n" %in% names(RUV)){
      k <- RUV$n
    } else {
      k <- 1
    }

    if ("empirical" %in% names(RUV)){
      ruv <- RUVSeq::RUVg(DESeq2::counts(dds, normalized = FALSE), RUV$empirical, k = k)

    } else if ("group" %in% names(RUV)){
      ruv <- RUVSeq::RUVs(DESeq2::counts(dds, normalized = FALSE), scIdx = RUVSeq::makeGroups(as.character(SummarizedExperiment::colData(dds)[[RUV$group]])), k = k)

    } else {
      stop("Error: Please provide either 'empirical' or 'group' to run RUVSeq!")
    }

    stopifnot(ncol(ruv$W) == k)
    if (!is.null(rownames(ruv$W))) stopifnot(all.equal(rownames(ruv$W), colnames(dds)))
    if (ncol(ruv$W) > 1){
      colnames(ruv$W) <- paste0("RUV", 1:ncol(ruv$W))
    } else {
      colnames(ruv$W) <- "RUV"
    }

    design_ruv <- data.frame(SummarizedExperiment::colData(dds) |> as.data.frame(), ruv$W)
    SummarizedExperiment::colData(dds) <- S4Vectors::DataFrame(design_ruv)
    fruv <- as.formula(paste0("~ . + ", paste0(colnames(ruv$W), collapse = " + ")))
    dds <- DESeq2::DESeqDataSet(dds, design = update.formula(old = formula, new = fruv))
  }

  # Explicit calculation of sizefactors ----
  if (!is.null(sizefactors)){
    dds$sizeFactor <- sizefactors
  } else if (!is.null(ctrlgenes)){
    if (!all(ctrlgenes %in% rownames(dds))){
      stop("Not all ctrlgenes present in data!")
    }
    dds <- DESeq2::estimateSizeFactors(dds, controlGenes = rownames(dds) %in% ctrlgenes)
  }


  # SVA ----
  if (length(SVA) > 0){
    if (!class(SVA) == "list") stop("Error: 'SVA' must be a named list.")
    if ("n" %in% names(SVA)){
      n <- SVA$n
    } else {
      n <- NULL
    }

    if (is.null(dds$sizeFactor)){
      dds <- DESeq2::estimateSizeFactors(dds)
    }

    if ("reduced" %in% names(SVA) & "formula" %in% class(SVA$reduced)){
      mm_full <- model.matrix(dds@design, SummarizedExperiment::colData(dds))
      mm_red <- model.matrix(SVA$reduced, SummarizedExperiment::colData(dds))
      svafit <- sva::svaseq(DESeq2::counts(dds, normalized = TRUE), mod = mm_full, mod0 = mm_red)

      if (svafit$n.sv > 0){
        if (!is.null(rownames(svafit$sv))) stopifnot(all.equal(rownames(svafit$sv), colnames(dds)))
        if (ncol(svafit$sv) > 1){
          colnames(svafit$sv) <- paste0("SV", 1:ncol(svafit$sv))
        } else {
          colnames(svafit$sv) <- "SV"
        }
        design_sva <- data.frame(SummarizedExperiment::colData(dds) |> as.data.frame(), svafit$sv)
        SummarizedExperiment::colData(dds) <- S4Vectors::DataFrame(design_sva)
        fsva <- as.formula(paste0("~ . + ", paste0(colnames(svafit$sv), collapse = " + ")))
        dds <- DESeq2::DESeqDataSet(dds, design = update.formula(old = formula, new = fsva))
      }

    } else {
      stop("Error: Please provide a reduced model formula to run SVA!")
    }

  }

  # Model fitting ----
  dds <- DESeq2::DESeq(dds, fitType = fitType, minReplicatesForReplace = minReplicatesForReplace, parallel = (ncores > 1), BPPARAM = bppar, ...)


  # Post-filtering ----
  if (!is.null(postfilter)){
    if (is.logical(postfilter)) stop("Error: Please provide gene IDs to 'postfilter'!")
    dds <- dds[intersect(postfilter, rownames(dds)),]
  }


  # Normalized counts ----
  results$normcounts <- DESeq2::counts(dds, normalized = TRUE)
  if (vst == TRUE) results$vst <- SummarizedExperiment::assays(DESeq2::vst(dds, blind = TRUE))[[1]]
  if (rlog == TRUE) results$rlog <- SummarizedExperiment::assays(DESeq2::rlog(dds, blind = TRUE))[[1]]


  # Contrasts ----
  if (!is.null(contrasts)){

    if (ihw == TRUE) filterFun <- IHW::ihw else filterFun <- NULL

    results$results <- lapply_fun(contrasts, BPPARAM = bppar, FUN = function(tmp){

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
      if (ordered == TRUE) mle <- dplyr::arrange(mle, pvalue)
      mle
    })


    # Make dataframes
    select_res <- function(contr, res, what){
      df <- res$results[[contr]]
      df$contrast <- contr
      df |> dplyr::select(gene, !!what, contrast)
    }

    results$log2FC <- lapply(names(contrasts), select_res, res = results, what = "log2FC") |> Reduce(f = rbind) |>
      tidyr::pivot_wider(names_from = "contrast", values_from = "log2FC") |> as.data.frame() |> col2rownames(gene)
    results$log2FC <-  results$log2FC[rownames(results$normcounts),, drop = FALSE]

    results$padj <- lapply(names(contrasts), select_res, res = results, what = "padj") |> Reduce(f = rbind) |>
      tidyr::pivot_wider(names_from = "contrast", values_from = "padj") |> as.data.frame() |> col2rownames(gene)
    results$padj <-  results$padj[rownames(results$normcounts),, drop = FALSE]
  }

  if (!is.null(lrt_reduced)){
    dds_rt <- DESeq2::DESeq(dds, test = "LRT", reduced = lrt_reduced,
                            fitType = fitType, minReplicatesForReplace = minReplicatesForReplace,
                            parallel = (ncores > 1), BPPARAM = bppar, ...)
    lrt <- DESeq2::results(dds)
    lrt$log2FoldChange <- NULL
    lrt$lfcSE <- NULL
    lrt$gene <- rownames(lrt)
    lrt <- as.data.frame(lrt)
    lrt <- dplyr::relocate(.data = lrt, gene)
    if (ordered == TRUE) lrt <- dplyr::arrange(lrt, pvalue)
    results$LRT <- lrt
  }

  # Results ----
  results$dds <- dds
  results$design <- as.data.frame(SummarizedExperiment::colData(dds))
  results

}










#' Differential expression analysis with LIMMA
#'
#' @param data Data
#' @param design Experimental design/colData
#' @param formula Formula
#' @param contrasts Named list of contrasts, specified as c(factor, level, reflevel) or c(numeric)
#' @param trend Trend
#' @param robust Robust
#' @param p.adj.method P.adj.method
#' @param do_log Log-transform the data
#' @param normalize Normalize
#' @param norm.method Normalization method
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' data <- getTestData()
#' runLIMMA(log2(data$data+1), design = data$design, formula = ~ group + cov1, contrasts = data$contrasts)
#' runLIMMA(data$data, do_log = TRUE, design = data$design, formula = ~ group + cov1, contrasts = list(cov1 = "cov1"))
runLIMMA <- function(data, design, formula = ~ 1, contrasts = NULL, trend = TRUE, robust = FALSE, p.adj.method = "fdr", do_log = FALSE, normalize = FALSE, norm.method = "vsn", ...){

  stopifnot(requireNamespace("limma", quietly = TRUE))

  if (ncol(data) != nrow(design)){
    stop("Error: The number of columns in data must match the number of rows in the design data frame.")
  }

  if (!is.null(colnames(data)) & !is.null(rownames(design))){
    if (!all(colnames(data) == rownames(design))){
      warning("Warning: Mismatch in 'colnames(data)' and 'rownames(design)'!")
    }
  }

  datamisc::colorcat("Limma differential expression analysis", col = "blue")
  datamisc::colorcat("-use 'do_log=TRUE' for log-transformation", col = "blue")
  datamisc::colorcat("-use 'normalize=TRUE' for normalization and 'norm.method' to set the method", col = "blue")

  if (do_log == TRUE){
    data <- log2(data + 1)
  }

  if (normalize == TRUE & tolower(norm.method) != "vsn"){
    if (do_log == TRUE){
      stop("Error: Do not combine log-transformation and vsn!")
    }
    data <- limma::normalizeVSN(data)
  } else if (normalize == TRUE){
    data <- limma::normalizeBetweenArrays(data, method = norm.method)
  }

  formula <- update(formula,  ~  0 + .)
  mm <- model.matrix(formula, design)
  if (qr(mm)$rank < ncol(mm)){
    warning("Model matrix is rank-deficient.")
  }

  contrasts_named <- lapply(contrasts, function(contr){
    if (length(contr) == 1){
      contr
    } else {
      paste(c(paste0(contr[1], contr[2]), paste0(contr[1], contr[3])), collapse = " - ")
    }
  })
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


  list(data = data, fit = efit, results = results)
}







#' Differential expression analysis with custom test
#'
#' @param data
#' @param design
#' @param contrasts
#' @param FUN
#' @param pull
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runTEST <- function(data, design, contrasts = NULL, FUN = t.test, pull = c(p = "p.value", stat = "statistic"), ...){

  res_na <- pull
  res_na[] <- NA
  results <- lapply(contrasts, function(contr){
    if (length(contr) > 1){
      contr_design <- design[design[[contr[1]]] %in% unlist(contr[-1]),,drop=FALSE]
    }
    contr_data <- data[,rownames(contr_design),drop=FALSE]
    group <- contr_design[[contr[1]]]

    contr_res <- apply(contr_data, 1, function(contr_data_row){
      g <- unique(group)
      x <- contr_data_row[naf(group == g[1])]
      y <- contr_data_row[naf(group == g[2])]
      if (length(unique(x[!is.na(x)])) >= 2 & length(unique(y[!is.na(y)])) >= 2){
        res <- FUN(x, y, ...)
        if (!is.null(pull)){
          sapply(pull, function(p) as.numeric(res[[p]]) )
        } else {
          res
        }
      } else {
        res_na
      }
    }) |> t() |> as.data.frame()
  })

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


  data[unique(na.omit(match(ids, rownames(data)))),, drop = FALSE] # bring into original order

}




#' Expand (demultiplex) rows with multiple ids
#'
#' @param data Data.frame
#' @param id_split List of ID vectors
#'
#' @return
#' @export
#'
#' @examples
#' x <- rmat()
#' rownames(x) <- paste0(rownames(x), ";s", 1:3)
#' expand(x, strsplit(rownames(x), split = ";"))
expand <- function(data, id_split){

  data <- as.data.frame(data)
  n <- sapply(id_split, length)
  data_expand <- data[n > 1,,drop=FALSE]
  id_expand <- id_split[n > 1]

  data_expand <- lapply(seq_along(id_expand), function(i){
    df <- data_expand[rep(i, length(id_expand[[i]])),,drop=FALSE]
    rownames(df) <- id_expand[[i]]
    df
  }) |> Reduce(f = dplyr::bind_rows)

  data <- dplyr::bind_rows(data[n == 1,,drop=FALSE], data_expand)
  if (!dplyr::setequal(rownames(data), unlist(id_split)) | any(duplicated(rownames(data)))){
    stop("Error in IDs!")
  }
  data[unlist(id_split),]
}







#' Get feature ranking vector from differential abundance testing results
#'
#' @description
#' Calculate ranking values to sort features from most upregulated to most downregulated for GSEA
#'
#' @param data
#' @param rank_by
#' @param rank_type
#'
#' @return
#' @export
#'
#' @examples
#' df <- data.frame(id = rgenes(n = 100), stat = rnorm(100), mean = runif(100))
#' df[1:50,-1] <- df[51:100,-1] # make duplicates
#' df |> get_feature_ranks(rank_by = c("stat", "mean"), rank_type = c("", ""), id_col="id")
#'
#' dds <- DESeq2::makeExampleDESeqDataSet(n = 10000, interceptMean = c(2,5))
#' res <- dds |> runDESeq2(formula = ~ condition, contrasts = list(BvsA = c("condition", "B", "A")))
#' res$results$BvsA |> get_feature_ranks(rank_by = c("pvalue", "baseMean"), rank_type = c("p", ""), direction_by = "log2FC")
get_feature_ranks <- function(data, rank_by = "stat", rank_type = "", id_col = NULL, direction_by = NULL, ...){

  datamisc::colorcat("Set 'rank_by' to define the main ranking statistic and secondary tie-breaking variables.", col = "blue")
  datamisc::colorcat("Use rank_type = '' to use 'rank_by' as ranking statistic", col = "blue")
  datamisc::colorcat("Use rank_type = 'p' to calculate the ranking statistic as -log10('rank_by') * sign('direction_by')", col = "blue")
  datamisc::colorcat("Use rank_type = '+' to weight 'rank_by' by sign('direction_by')", col = "blue")

  data <- as.data.frame(data)
  stopifnot(length(rank_by) == length(rank_type))
  keep <- unique(c(rank_by, id_col, direction_by))
  if (!all(keep %in% colnames(data))){
    warning(paste0("Column ", setdiff(keep, colnames(data)), " not found in data."))
  }
  keep <- keep[keep %in% colnames(data)]
  df <- data[,keep,drop=FALSE]
  df <- df[!is.na(df[[1]]),,drop=FALSE]

  ## Transform values so that they can be sorted from highest to smallest

  # check if negative values are present

  for (i in seq_along(rank_by)){

    # "+" - sort from highest to smallest, needs direction indicator
    if (rank_type[i] == "+"){
      if (any(df[[rank_by[i]]] < 0)){ stop(paste0("Error: Encountered negative values in column ", rank_by[i])) }
      df[[rank_by[i]]] <- df[[rank_by[i]]] * sign(df[[direction_by]])

    # "p" - sort from smallest to highest, needs direction indicator
    } else if (rank_type[i] == "p"){
      if (any(df[[rank_by[i]]] < 0)){ stop(paste0("Error: Encountered negative values in column ", rank_by[i])) }
      df[[rank_by[i]]] <- -log10(pmax(df[[rank_by[i]]], .Machine$double.xmin)) * sign(df[[direction_by]])
    }
  }

  if (!is.null(id_col)){
    df[[id_col]] <- factor(df[[id_col]], levels = sort(unique(df[[id_col]]), decreasing = TRUE), ordered = TRUE)
  }

  df <- dplyr::arrange_all(df)
  df <- df[rev(1:nrow(df)),,drop = FALSE]

  ## Subranking
  df$rankid <- as.character(df[[1]])
  df$rankid[is.na(df$rankid)] <- "na"
  df$dup <- df$rankid %in% df$rankid[duplicated(df$rankid)]
  df$newstat <- df[[1]]
  if (!is.null(id_col)){
    df[[id_col]] <- as.character(df[[id_col]])
  }

  for (rankid in unique(subset(df, dup == TRUE)$rankid)){

    tmp <- df[df$rankid == rankid,,drop = FALSE]
    ix <- range(which(rankid == df$rankid))

    tmp2 <- df[c(ix[1]-1, ix[2]+1),,drop = FALSE]
    tmp2[is.na(tmp2)] <- mean(abs(diff(df[[1]])))
    maxdiff <- min(abs(tmp2[[1]] - mean(tmp[[1]]))) * 0.4

    tmpranks <- rev(1:nrow(tmp))
    tmpranks <- (tmpranks - median(tmpranks)) / max(abs(tmpranks))

    tmp$newstat <- tmp[[1]] + maxdiff * tmpranks
    df[df$rankid == rankid,]$newstat <- tmp$newstat
  }

  if (is.null(id_col)){
    ranks <- setNames(df$newstat, rownames(df))
  } else {
    ranks <- setNames(df$newstat, df[[id_col]])
  }

  ranks
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
#' @param rank_type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' data.frame(id = rgenes(n = 100), stat = rnorm(100)) |> col2rownames() |> runGSEA(genesets = getMSigDB(collections = "H"))
#'
#' dds <- DESeq2::makeExampleDESeqDataSet(n = 10000, interceptMean = c(2,5))
#' rownames(dds) <- rgenes(n = 10000)
#' res <- dds |> runDESeq2(formula = ~ condition, contrasts = list(BvsA = c("condition", "B", "A")))
#' res$results$BvsA |> runfGSEA(genesets = getMSigDB(collections = "H"))
runGSEA <- function(data, genesets = NULL, rank_by = "stat", rank_type = "", id_col = NULL, direction_by = NULL, min_size = 3, max_size = 2000, seed = 0, as.df = TRUE, ...){

  stopifnot(requireNamespace("clusterProfiler"))
  stopifnot(requireNamespace("fgsea"))

  if (is.null(genesets)) {
    genesets <- NA
  }
  if (!"data.frame" %in% class(genesets)) {
    genesets <- datamisc::convertGeneSets(genesets, to = "data.frame")
  }

  ranks <- datamisc::get_feature_ranks(data, rank_by = rank_by, rank_type = rank_type, id_col = id_col, direction_by = direction_by)

  set.seed(seed)
  results <- clusterProfiler::GSEA(ranks, seed = TRUE, eps = 0, minGSSize = min_size, maxGSSize = max_size, TERM2GENE = genesets, pvalueCutoff = 1, pAdjustMethod = "fdr", ...)

  if (as.df == TRUE) {
    results <- as.data.frame(results)
    results <- dplyr::select(results, -Description)
    results <- dplyr::mutate(results, core_enrichment = gsub("/", ", ", core_enrichment))
    results <- dplyr::rename(results, term = ID, ES = enrichmentScore, padj = p.adjust, qval = qvalue, gene_set_size = setSize, enriched_genes = core_enrichment)
    results$n_enriched_genes <- sapply(strsplit(results$enriched_genes, split = ", ", fixed = TRUE), FUN = length)
    results <- dplyr::relocate(results, n_enriched_genes, .after = gene_set_size)
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
#' data.frame(id = rgenes(n = 100), stat = rnorm(100)) |> col2rownames() |> runfGSEA(genesets = getMSigDB(collections = "H"))
#'
#' dds <- DESeq2::makeExampleDESeqDataSet(n = 10000, interceptMean = c(2,5))
#' rownames(dds) <- rgenes(n = 10000)
#' res <- dds |> runDESeq2(formula = ~ condition, contrasts = list(BvsA = c("condition", "B", "A")))
#' res$results$BvsA |> runfGSEA(genesets = getMSigDB(collections = "H"))
runfGSEA <- function(data, genesets, rank_by = "stat", rank_type = "", id_col = NULL, direction_by = NULL, min_size = 3, max_size = 2000, seed = 123){

  stopifnot(requireNamespace("fgsea"))

  if ("data.frame" %in% class(genesets)){
    genesets <- convertGeneSets(genesets, to = "list")
  }

  data <- as.data.frame(data)
  ranks <- datamisc::get_feature_ranks(data, rank_by = rank_by, rank_type = rank_type, id_col = id_col, direction_by = direction_by)

  set.seed(seed)
  results <- fgsea::fgsea(genesets, sort(ranks), minSize = min_size, maxSize = max_size, eps = 0)
  results <- as.data.frame(results)
  results <- dplyr::rename(results, term = pathway)
  results <- dplyr::arrange(results, pval)
  results
}



#' Single-sample gene set enrichment with GSVA
#'
#' @param data Expression data
#' @param genesets Genesets
#' @param method
#' @param kcdf
#' @param min_size
#' @param max_size
#' @param assay
#' @param annotation
#' @param params
#' @param ncores
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' data <- getTestData()$data
#' genesets = getMSigDB("H")
#' rownames(data) <- sample(unique(unlist(genesets)), nrow(data))
#' runGSVA(data, genesets)
runGSVA <- function(data, genesets, method = "gsva", kcdf = "Gaussian", min_size = 1, max_size = Inf, assay = NA_character_, annotation = NULL, params = NULL, ncores = 1, ...){

  stopifnot(requireNamespace("GSVA"))

  if ("data.frame" %in% class(genesets)){
    genesets <- convertGeneSets(genesets, to = "list")
  }

  ncores <- min(ncores, parallel::detectCores())

  # method params
  method <- tolower(method)
  if (method == "gsva" & is.null(params)){
    params <- GSVA::gsvaParam(exprData = data, geneSets = genesets, minSize = min_size, maxSize = max_size, assay = assay, annotation = annotation, kcdf = kcdf)
  }
  if (method == "ssgsea" & is.null(params)){
    params <- GSVA::ssgseaParam(exprData = data, geneSets = genesets, minSize = min_size, maxSize = max_size, assay = assay, annotation = annotation)
  }
  if (method == "plage" & is.null(params)){
    params <- GSVA::plageParam(exprData = data, geneSets = genesets, minSize = min_size, maxSize = max_size, assay = assay, annotation = annotation)
  }
  if (method == "zscore" & is.null(params)){
    params <- GSVA::zscoreParam(exprData = data, geneSets = genesets, minSize = min_size, maxSize = max_size, assay = assay, annotation = annotation)
  }

  # run method
  res <- GSVA::gsva(param = params, ...)
  as.data.frame(res)
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









