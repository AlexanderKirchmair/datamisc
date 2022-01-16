



#' Ggplot2 barplot
#'
#' @param x (named) numeric vector
#' @param color single color or vector of colors
#' @param ... arguments passed to geom_bar()
#'
#' @export
#'
#' @examples
#' ggbar(1:10)
#' ggbar(setNames(1:5, LETTERS[1:5]), color = rgb((1:5)/5, 0.4, 0.4), width = 0.5)
ggbar <- function(x, color = NULL, ...){

  df <- data.frame(row.names = names(x), x = seq(x), y = x)
  if (!is.null(names(x))) df$x <- factor(names(x), ordered = TRUE, levels = names(x))

  fill <- NULL
  if (length(color) > 1) df$fill <- color

  gg <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x, y, fill = fill)) +
    theme_basic() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(0)) +
    ggplot2::xlab("") +
    ggplot2::ylab("")

  if (length(color) == 1){
    gg <- gg + ggplot2::geom_bar(stat = "identity", fill = color, ...)
  } else {
    gg <- gg + ggplot2::geom_bar(stat = "identity", ...)
  }

  if (length(color) > 1) gg <- gg + ggplot2::scale_fill_identity(guide = "none")

  gg
}








#' Plot experimental design matrix
#'
#' @param design
#' @param columns
#' @param colors
#' @param label
#' @param legend
#' @param fontsize
#' @param nacol
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ggdesign <- function(design, columns = NULL, colors = NULL, label = NULL, legend = NULL, fontsize = 15, nacol = "grey70", ...){

  cols <- rlang::enquo(columns)

  if (is.null(label)) label <- nrow(design) <= 30
  if (is.null(legend)) legend <- nrow(design) > 30

  if (!rlang::quo_is_null(cols)) design <- dplyr::select(design, !!cols)

  tmp <- design[,!colnames(design) %in% names(colors), drop = FALSE]
  colors <- c(colors, getColors(tmp))

  design$Sample <- rownames(design)
  df <- tidyr::pivot_longer(design, cols = -Sample, names_to = "Factor")
  df$Factor <- factor(df$Factor, ordered = TRUE, levels = colnames(design))
  df$Sample <- factor(df$Sample, ordered = TRUE, levels = rev(rownames(design)))

  gg <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x = Factor, y = Sample, fill = value, label = value, ...)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(colour = "black", size = fontsize)) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::xlab("") + ggplot2::ylab("")

    if (!is.null(colors)){
      gg <- gg + ggplot2::scale_fill_manual(values = unlist(setNames(colors[names(colors) %in% colnames(design)], NULL), recursive = FALSE), na.value = nacol, guide = "none")
    }


  if (label == TRUE) gg <- gg + ggplot2::geom_text(size = fontsize/ggplot2::.pt, alpha = 0.6)
  if (legend == FALSE) gg <- gg + ggplot2::theme(legend.position = "none")

  gg
}









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
ggpca <- function(data, design = NULL, mapping =  ggplot2::aes(), center = TRUE, scale = TRUE, na = NULL, n = NULL, label = FALSE, digits = 2, colors = NULL, size = 3.5, return_data = FALSE, ...){

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
    ggplot2::theme(panel.border =  ggplot2::element_rect(color = "black", fill = NA, size = 1),
                   axis.text =  ggplot2::element_text(color = "black"),
                   axis.ticks =  ggplot2::element_line(color = "black")) +
    ggplot2::geom_point(size = size) +
    ggplot2::xlab(paste0("PC1 (", round(100*variance_explained[1], digits), "%)")) +
    ggplot2::ylab(paste0("PC2 (", round(100*variance_explained[2], digits), "%)")) +
    ggplot2::scale_x_continuous(expand =  ggplot2::expansion(mult = c(0.2,0.2))) +
    ggplot2::scale_y_continuous(expand =  ggplot2::expansion(mult = c(0.2,0.2)))

  if (label == TRUE) gg <- gg +  ggplot2::geom_text_repel(size = 5, show.legend = F, min.segment.length = 2)
  if (!is.null(colors)) gg <- gg +  ggplot2::scale_color_manual(values = colors[[rlang::as_name(base_aes[["colour"]])]])

  gg

}








#' Save ggplot as pdf
#'
#' @param gg
#' @param file
#' @param width
#' @param height
#' @param dpi
#' @param units
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ggpdf <- function(gg, file, width = 3000, height = 2500, dpi = 300, units = "px", ...){
  if (length(file) > 1) file <- do.call(file.path, as.list(file))
  if (nat(baseext(file) != "pdf")) file <- paste0(file, ".pdf")
  ggplot2::ggsave(filename = file, plot = gg, width = width, height = height, dpi = 300, units = units, device = "pdf", ...)
  invisible(gg)
}


#' Save ggplot as png
#'
#' @param gg
#' @param file
#' @param width
#' @param height
#' @param dpi
#' @param units
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ggpng <- function(gg, file, width = 3000, height = 2500, dpi = 300, units = "px", ...){
  if (length(file) > 1) file <- do.call(file.path, as.list(file))
  if (nat(baseext(file) != "png")) file <- paste0(file, ".png")
  ggplot2::ggsave(filename = file, plot = gg, width = width, height = height, dpi = 300, units = units, device = "png", type = "cairo", ...)
  invisible(gg)
}








gg_getLimits <- function(gg){
  gb <- ggplot2::ggplot_build(gg)
  x <- gb$layout$panel_params[[1]]$x.range
  y <- gb$layout$panel_params[[1]]$y.range
  list(x = x, y = y)
}


gg_getGeoms <- function(gg){
  sapply(gg$layers, function(tmp) class(get("geom", tmp))[1] )
}



#' Ggplot with errorbars
#'
#' @param data
#' @param mapping
#' @param type
#' @param barwidth
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ggstat <- function(data, mapping, type = sd, barwidth = 0.1, ...){

  # input
  data <- data.frame(data)
  x <- mapping[["x"]]
  y <- mapping[["y"]]
  group_by <- mapping[["colour"]]
  type <- rlang::enquo(type) # must be sd, var or se

  # summarize data
  statsdf <- subset(data, !is.na(rlang::as_name(x)) & !is.na(rlang::as_name(y))) %>%
    group_by(!!group_by, !!x) %>% summarize(mean = mean(!!y, na.rm = TRUE),
                                            var = var(!!y, na.rm = TRUE),
                                            sd = sd(!!y, na.rm = TRUE),
                                            n = n())

  statsdf %<>% mutate(se = sqrt(var/n) )
  statsdf %<>% mutate(ymin = mean - !!type, ymax = mean + !!type)

  # ggplot
  ggs <- ggplot2::ggplot(data = statsdf, mapping = aes(x = !!x, y = mean, colour = !!group_by)) + # mean values as default data
    ggplot2::geom_point(data = data, mapping = mapping) +
    geom_errorbar(aes(ymin = ymin, max = ymax), width = barwidth, ...) +
    ylab(rlang::as_name(y))

  getbw <- function(ggs, barwidth){diff(gg_getLimits(ggs)$x) * barwidth}
  ggs$layers[[which(gg_getGeoms(ggs) == "GeomErrorbar")]]$geom_params$width <- getbw(ggs, barwidth) # consistent bar widths
  ggs
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
#' cxheatmap(rmat(50, 30), coldf = data.frame(row.names = colnames(rmat(50, 30)), group = sample(LETTERS[1:5], size = 30, replace = TRUE)))
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
  heatdata <- data.matrix(heatdata)
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
  addcol <- NULL
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
      border <- grid::gpar(col = rgb(1,1,1), lwd = grid::unit(1, "pt"))
    } else {
      border <- grid::gpar(col = NA)
    }
  } else {
    if (any(border == TRUE)){
      border <- grid::gpar(col = rgb(1,1,1), lwd = grid::unit(1, "pt"))
    } else if (length(border) > 1){
       border <- grid::gpar(col = border[is.na(as.numeric(border))], lwd = grid::unit(as.numeric(border)[!is.na(as.numeric(border))], "pt"))
    } else {
       border <- grid::gpar(col = NA)
    }
  }




  ### Annotations ----

  # Legends
  # see ?`color_mapping_legend,ColorMapping-method`
  legend_params <- list(title_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
                        legend_height = grid::unit(0.2, "npc"),
                        border = legend_border,
                        labels_gp = grid::gpar(fontsize = fontsize))


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
      if (naf(cellmat[i,j] == TRUE)){ grid::grid.points(x, y, pch = markshape, size = unit(marksize, "pt")) }
    }
  }



  ### Heatmap ----

  dimnames(heatdata) <- lapply(dimnames(heatdata), function(x) cutstr(x, maxchar = maxchar))

  if (is.null(rowcex) & nrow(heatdata) > 10*ncol(heatdata)) rowcex <- 1/log10(nrow(heatdata))
  if (is.null(rowcex)) rowcex <- 1
  if (is.null(colcex)) colcex <- 1

  hm <- ComplexHeatmap::Heatmap(name = title,
                                matrix = heatdata,
                                row_names_max_width = grid::unit(0.3, "npc"),
                                column_names_max_height = grid::unit(0.3, "npc"),
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

  args <- list(df = df,
               name = side,
               annotation_label = colnames(df),
               which = ifelse(side %in% c("top", "bottom"), "column", "row"),
               col = colors[names(colors) %in% colnames(df)],
               show_annotation_name = TRUE,
               gap = grid::unit(gap, "cm"),
               border = FALSE,
               gp = grid::gpar(col = rgb(0,0,0)),
               annotation_name_gp = grid::gpar(fontsize = fontsize, fontface = "bold"),
               simple_anno_size_adjust = TRUE,
               annotation_name_side = ifelse(side %in% c("top", "bottom"), "right", "top"),
               annotation_legend_param = legend_params,
               show_legend = legend,
               ...)

  do.call(ComplexHeatmap::HeatmapAnnotation, args)

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








