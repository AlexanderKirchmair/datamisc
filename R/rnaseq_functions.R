


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





#' Import and plot preprocessing statistics from nf-core/rnaseq pipeline
#'
#' @param multiqc_dir
#' @param design
#' @param ignore
#'
#' @return
#' @export
#'
#' @examples
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




