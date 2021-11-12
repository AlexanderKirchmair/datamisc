



#' Principal component analysis and plotting of data
#'
#' @description
#' Quick plotting of PCA results.
#' For advanced analysis, see e.g. browseVignettes("PCAtools").
#'
#' @param data
#' @param design
#' @param mapping
#' @param center
#' @param scale
#' @param na
#' @param n
#' @param label
#' @param digits
#' @param colors
#' @param size
#' @param return_data
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ggpca <- function(data, design = NULL, mapping = aes(), center = TRUE, scale = TRUE, na = NULL, n = NULL, label = FALSE, digits = 2, colors = NULL, size = 3.5, return_data = FALSE, ...){

  # tidy data case


  # annotation data
  if (is.null(design)) design <- data.frame(row.names = rownames(data))

  # missing values
  if (any(is.na(data))){
    if (na == "impute"){
      ncp <-  missMDA::estim_ncpPCA(data, ncp.min = 2, ncp.max = max(round(ncol(data)/3), 3))
      data <- missMDA::imputePCA(data, ncp = ncp$ncp)$completeObs
    } else if (na == "omit"){
      data <- na.omit(data)
    } else if (is.numeric(na)){
      data[is.na(data)] <- na
    } else {
      stop("Please select a method for dealing with missing values!")
    }
  }


  # zero-variance data
  vars <- apply(data, 1, var, na.rm = TRUE)
  data <- data[abs(vars) > .Machine$double.eps & !is.na(vars),]

  # top-variance features
  if (!is.null(n)){
    data <- data[,order(apply(data, 2, var), decreasing = TRUE)[1:n]]
  }

  # run pca
  pcares <- stats::prcomp(t(data), center = center, scale. = scale, ...)
  if (return_data == TRUE) return(pcares)

  # plotting data
  pca_df <- data.frame(pcares$x, design[rownames(pcares$x),,drop = FALSE])
  variance_explained <- summary(pcares)$importance["Proportion of Variance", ]
  pca_df$label <- rownames(pca_df)

  # ggplot
  base_aes <- ggplot2::aes(x = PC1, y = PC2, label = label)
  base_aes[names(mapping)] <- mapping

  gg <- ggplot2::ggplot(data = pca_df, mapping = base_aes) +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
                   axis.text = element_text(color = "black"),
                   axis.ticks = element_line(color = "black")) +
    ggplot2::geom_point(size = size) +
    ggplot2::xlab(paste0("PC1 (", round(100*variance_explained[1], digits), "%)")) +
    ggplot2::ylab(paste0("PC2 (", round(100*variance_explained[2], digits), "%)")) +
    ggplot2::scale_x_continuous(expand = expansion(mult = c(0.2,0.2))) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.2,0.2)))

  if (label == TRUE) gg <- gg + geom_text_repel(size = 5, show.legend = F, min.segment.length = 2)
  if (!is.null(colors)) gg <- gg + scale_color_manual(values = colors[[rlang::as_name(base_aes[["colour"]])]])

  gg

}









### Heatmap functions ----





#' Heatmap plotting with ComplexHeatmap
#'
#' @param data matrix or dataframe
#' @param rowdf dataframe with row annotations
#' @param coldf dataframe with column annotations
#' @param scale z-scale 'rows' or 'cols'
#' @param cluster_rows cluster rows
#' @param cluster_cols cluster columns
#' @param rowdf_side left/right
#' @param coldf_side top/bottom
#' @param rowdf_legend show rowdf legend
#' @param coldf_legend show coldf legend
#' @param fontsize base fontsize
#' @param rowcex sizefactor for row annotations
#' @param colcex sizefactor for col annotations
#' @param heatpal color palette ('circlize::colorRamp2')
#' @param border cell border
#' @param title title of colorscale legend
#' @param colors list of colors for annotations
#' @param mat logical matrix indicating where to draw marks
#' @param markshape shape of cell marks
#' @param marksize size of cell marks
#' @param na_col color of NA values
#' @param maxchar max. length of colnames/rownames
#' @param legend_border border of legends
#' @param inf handling of Inf values for clustering
#' @param na handling of NA values for clustering
#' @param markoob mark out-of-bounds values
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' cxheatmap(rmat(50, 30), coldf = data.frame(row.names = colnames(rmat(50, 30)), group = sample(LETTERS[1:5], size = 30, replace = T)))
#' cxheatmap(rmat(50, 30), border = c("black", 2))
#'
#' @seealso https://jokergoo.github.io/ComplexHeatmap-reference/book/
#'
cxheatmap <- function(data, rowdf = NULL, coldf = NULL, scale = FALSE, cluster_rows = NULL, cluster_cols = NULL,
                      rowdf_side = "left", coldf_side = "top", rowdf_legend = TRUE, coldf_legend = TRUE,
                      legend_border = "black",
                      fontsize = 12, rowcex = NULL, colcex = 1,
                      heatpal = NULL, border = NULL, title = NULL, colors = NULL,
                      inf = F, na = 0, mat = NULL, markoob = FALSE, markshape = 4, marksize = NULL, na_col = "grey", maxchar = 35, ...){


  ### Data ----
  datacall <- substitute(data)
  if (is.null(title)){
    if (grepl("row|col", scale, ignore.case = TRUE)){
      title <- "z-score"
    } else {
      title <- deparse1(datacall)
    }
  }
  if (title == FALSE) title <- " "

  heatdata <- eval(datacall, envir = parent.frame())
  heatdata <- matScale(heatdata, rows = grepl("row", scale, ignore.case = TRUE), cols = grepl("col", scale, ignore.case = TRUE))

  ### Clustering ----
  clust <- clusterData(heatdata, rows = cluster_rows, cols = cluster_cols, inf = inf, na = na)
  if (is.null(clust$rows)) clust$rows <- FALSE
  if (is.null(clust$cols)) clust$cols <- FALSE


  ### Colors ----

  # scale colors
  if (is.null(heatpal)){
    heatpal_colors <- getColorScale(heatdata)
    heatpal <- circlize::colorRamp2(breaks = heatpal_colors, colors = names(heatpal_colors))
  }

  # annotation colors
  docol <- setdiff(unlist(lapply(list(coldf, rowdf), colnames)), names(colors))
  if (length(docol) > 0){
    if (!is.null(coldf)) all <- coldf
    if (!is.null(rowdf)) all <- rowdf
    if (!is.null(coldf) & !is.null(rowdf)) all <- dplyr::full_join(coldf, rowdf, by = character())
    addcol <- getColors(all[,docol,drop = FALSE])
  }

  colors <- c(colors, addcol)
  colors <- lapply(colors, function(tmp) tmp[!is.na(tmp) & !is.na(names(tmp))])

  # cell border
  if (is.null(border)){
    if (nrow(heatdata) < 100 & ncol(heatdata) < 100){
      border <- grid::gpar(col = rgb(1,1,1), lwd = unit(1, "pt"))
    } else {
      border <- grid::gpar(col = NA)
    }
  } else {
    if (any(border == TRUE)){
      border <- grid::gpar(col = rgb(1,1,1), lwd = unit(1, "pt"))
    } else if (length(border) > 1){
       border <- grid::gpar(col = border[is.na(as.numeric(border))], lwd = unit(as.numeric(border)[!is.na(as.numeric(border))], "pt"))
    } else {
       border <- grid::gpar(col = NA)
    }
  }




  ### Annotations ----

  # Legends
  # see ?`color_mapping_legend,ColorMapping-method`
  legend_params <- list(title_gp = gpar(fontsize = fontsize, fontface = "bold"),
                        legend_height = unit(0.2, "npc"),
                        border = legend_border,
                        labels_gp = gpar(fontsize = fontsize))


  # Row annotation
  rowAnn <- NULL
  if (!is.null(rowdf)) rowAnn <- getCXanno(df = rowdf[rownames(heatdata),, drop = FALSE],
                                           colors = colors,
                                           side = rowdf_side,
                                           legend = rowdf_legend,
                                           legend_params = legend_params)

  # Column annotation
  colAnn <- NULL
  if (!is.null(coldf)) colAnn <- getCXanno(coldf[colnames(heatdata),, drop = FALSE],
                                           colors = colors,
                                           side = coldf_side,
                                           legend = coldf_legend,
                                           legend_params = legend_params)


  # Cell annotation

  if (markoob == TRUE & is.null(mat)){
    mat <- matrix(data = FALSE, nrow = nrow(heatdata), ncol = ncol(heatdata), dimnames = dimnames(heatdata))
    mat[heatdata < min(heatpal_colors)] <- TRUE
    mat[heatdata > max(heatpal_colors)] <- TRUE
  }

  if (is.null(marksize)) marksize <- fontsize * 0.6
  cellFUN <- NULL
  if (!is.null(mat)){
    if (!is.logical(mat)) stop("'Mat' must be a logical indicator of whether cells should be marked!")
    cellmat <- mat[rownames(heatdata), colnames(heatdata)]
    cellFUN <- function(j, i, x, y, width, height, fill){
      if (naf(cellmat[i,j] == TRUE)){ grid.points(x, y, pch = markshape, size = unit(marksize, "pt")) }
    }
  }



  ### Heatmap ----

  dimnames(heatdata) <- lapply(dimnames(heatdata), function(x) cutstr(x, maxchar = maxchar))

  if (is.null(rowcex) & nrow(heatdata) > 10*ncol(heatdata)) rowcex <- 1/log10(nrow(heatdata))
  if (is.null(rowcex)) rowcex <- 1
  if (is.null(colcex)) colcex <- 1

  hm <- ComplexHeatmap::Heatmap(name = title,
                                matrix = heatdata,
                                row_names_max_width = unit(0.3, "npc"),
                                column_names_max_height = unit(0.3, "npc"),
                                column_title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
                                rect_gp = border,
                                na_col = na_col,
                                left_annotation = rowAnn,
                                top_annotation = colAnn,
                                row_names_gp = grid::gpar(fontsize = fontsize * rowcex),
                                column_names_gp = grid::gpar(fontsize = fontsize * colcex),
                                col = heatpal,
                                heatmap_legend_param = legend_params,
                                cluster_rows = clust$rows,
                                cluster_columns = clust$cols,
                                cell_fun = cellFUN,
                                ...)

  hm
}




getCXanno <- function(df = NULL, side = "top", colors = NULL, fontsize = 12, gap = 0, legend = FALSE, legend_params = NULL, ...){

  if (is.null(df)) return(list())

  ComplexHeatmap::HeatmapAnnotation(df = df,
                                    name = side,
                                    annotation_label = colnames(df),
                                    which = ifelse(side %in% c("top", "bottom"), "column", "row"),
                                    col = colors[names(colors) %in% colnames(df)],
                                    show_annotation_name = TRUE,
                                    gap = unit(gap, "cm"),
                                    border = FALSE,
                                    gp = grid::gpar(col = rgb(0,0,0)),
                                    annotation_name_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
                                    simple_anno_size_adjust = TRUE,
                                    annotation_name_side = ifelse(side %in% c("top", "bottom"), "right", "top"),
                                    annotation_legend_param = legend_params,
                                    show_legend = legend,
                                    ...)

}






getColorScale <- function(data, ...){

  data <- data.matrix(data)
  data[!is.finite(data)] <- NA
  lims <- range(data, na.rm = TRUE)

  if (min(data, na.rm = TRUE) >= 0){ # positive only
    col <- setNames(c(0, lims[2]), c("white", "red"))

  } else {
    col <- setNames(c(lims[1], 0, lims[2]), c("blue", "white", "red"))

  }

  col
}








