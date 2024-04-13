




#' Barplots with ggplot2
#'
#' @param data
#' @param aes
#' @param colors
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' ggbar(setNames(1:5, LETTERS[1:5]))
#' ggbar(mtcars$mpg)
ggbar <- function(data, aes = NULL, colors = NULL, x.cex = 1, rot = NULL, ...){

  if (is.null(dim(data))){
    # vector input
    df <- data.frame(row.names = names(data), x = seq_along(data), y = data)
    if (!is.null(names(data))) df$x <- factor(names(data), ordered = TRUE, levels = names(data))
    if (is.null(colors)) colors <- "#76a3b8"
    aes <- ggplot2::aes(x, y)

  } else {
    # data.frame input
    df <- data
    if (is.null(df$x)) df$x <- 1:nrow(df)
    if (is.null(df$y)) df$y <- df[[2]]
    if (is.null(aes)) aes <- ggplot2::aes(x, y)
  }

  # plotting
  max_chars <- max(nchar(as.character(df$x)))
  if (is.null(rot)){
    rot <- max_chars > 3
  }


  gg <- ggplot2::ggplot(data = df, mapping = aes) +
    ggplot2::`%+replace%`(theme_basic(),
                          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = ifelse(rot, 90, 0), size = ggplot2::rel(x.cex),
                                                                             vjust = ifelse(rot, 0.5, 0.5),
                                                                             hjust = ifelse(rot, 1, 0.5))) ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(0)) +
    ggplot2::xlab("") +
    ggplot2::ylab("")

  if (length(colors) == 1){
    gg <- gg + ggplot2::geom_bar(stat = "identity", fill = colors, ...)
  } else {
    gg <- gg + ggplot2::geom_bar(stat = "identity", ...)
  }

  if (length(colors) > 1) gg <- gg + ggplot2::scale_fill_manual(values = colors)

  gg
}






#' Plot colData/experimental design matrix
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

  # select top-variance features
  if (!is.null(n)){
    if (n <= 1) n <- round(nrow(data) * n)
    vars <- apply(data, 1, var)
    data <- data[order(abs(vars), decreasing = TRUE)[1:n],]
  }

  # zero-variance data
  vars <- apply(data, 1, var, na.rm = TRUE)
  data <- data[abs(vars) > .Machine$double.eps & !is.na(vars),]

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

  if (label == TRUE) gg <- gg + ggrepel::geom_text_repel(size = 5, show.legend = F, min.segment.length = 2, seed = 123)
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









#' Save plot to file
#'
#' @param p plot
#' @param file filename
#' @param dev device
#' @param width
#' @param height
#' @param dpi
#' @param units
#' @param bg
#' @param gg
#' @param PRINTFUN
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
saveplot <- function(p, file = NULL, dev = "png", width = 3000, height = 2500, dpi = 300, units = "px", bg = "white", ggsave = TRUE, PRINTFUN = print, ...){


  stopifnot(dev %in% c("png", "pdf", "tiff", "svg", "jpeg", "jpg", "cairo_ps", "ps"))

  # filename
  if (is.null(file)) file <- rlang::as_name(rlang::enquo(p))
  if (length(file) > 1){ file <- do.call(file.path, as.list(file)) }
  if (nat(baseext(file) != dev)) file <- paste0(file, ".", dev)



  get_inches <- function(x, dpi, units){
    if (units == "px"){
      x <- x/dpi
    }
    if (units == "cm"){
      x <- x/dpi / 2.54
    }
    x
  }



  # GGSAVE WRAPPER
  # only for single ggplots
  if ("gg" %in% class(p) & ggsave == TRUE & !(dev %in% c("pdf", "ps", "cairo_ps")) ){
    ggplot2::ggsave(filename = file,
                    plot = p,
                    device = dev,
                    width = width,
                    height = height,
                    units = units,
                    dpi = dpi,
                    bg = bg,
                    ...)

  } else {

    # PNG
    if (dev == "png"){
      if (is.null(bg)) bg <- "white"
        grDevices::png(filename = file,
                       width = width,
                       height = height,
                       units = units,
                       bg = bg,
                       res = dpi,
                       type = "cairo",
                       ...)
          PRINTFUN(p)
        dev.off()
    }


    # PDF
    if (dev == "pdf"){
      # multi-page pdf if list of plots (with lenght > 1)
      if (!"list" %in% class(p)) p <- list(p)
      grDevices::pdf(file = file,
                     width = get_inches(width, dpi = dpi, units = units),
                     height = get_inches(height, dpi = dpi, units = units),
                     bg = bg,
                     onefile = TRUE,
                     ...)
        tmp <- lapply(p, PRINTFUN)
      dev.off()
    }


    # SVG
    if (dev == "svg"){
      grDevices::svg(filename = file,
                     width = get_inches(width, dpi = dpi, units = units),
                     height = get_inches(height, dpi = dpi, units = units),
                     onefile = TRUE,
                     ...)
        PRINTFUN(p)
      dev.off()
    }

    # TIFF
    if (dev == "tiff"){
      grDevices::tiff(filename = file,
                      width = width,
                      height = height,
                      res = dpi,
                      bg = bg,
                      type = "cairo",
                      units = units,
                      ...)
        PRINTFUN(p)
      dev.off()
    }

    # JPG
    if (dev %in% c("jpeg", "jpg")){
      grDevices::jpeg(filename = file,
                      width = width,
                      height = height,
                      res = dpi,
                      bg = bg,
                      type = "cairo",
                      units = units,
                      ...)
        PRINTFUN(p)
      dev.off()
    }

    # POSTSCRIPT
    if (dev %in% c("ps", "cairo_ps")){
      grDevices::cairo_ps(filename = file,
                          width = get_inches(width, dpi = dpi, units = units),
                          height = get_inches(height, dpi = dpi, units = units),
                          fallback_resolution = dpi,
                          bg = bg,
                          ...)
        PRINTFUN(p)
      dev.off()
    }

  }


  invisible(p)
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
                      legend_border = "black", anno_border = "black",
                      fontsize = 12, rowcex = NULL, colcex = 1,
                      rownames_width = 0.3, colnames_width = 0.3,
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
                                           anno_border = anno_border,
                                           side = rowdf_side,
                                           legend = rowdf_legend,
                                           legend_params = legend_params)

  # Column annotation
  colAnn <- NULL
  if (!is.null(coldf)) colAnn <- getCXanno(coldf[colnames(heatdata),, drop = FALSE],
                                           colors = colors,
                                           anno_border = anno_border,
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
                                row_names_max_width = grid::unit(rownames_width, "npc"),
                                column_names_max_height = grid::unit(colnames_width, "npc"),
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




getCXanno <- function(df = NULL, side = "top", colors = NULL, fontsize = 12, gap = 0, legend = FALSE, anno_border = NULL, legend_params = NULL, ...){

  if (is.null(df)) return(list())

  if (is.null(anno_border)){
    gp <- grid::gpar()
  } else {
    gp <- grid::gpar(col = anno_border)
  }

  args <- list(df = df,
               name = side,
               annotation_label = colnames(df),
               which = ifelse(side %in% c("top", "bottom"), "column", "row"),
               col = colors[names(colors) %in% colnames(df)],
               show_annotation_name = TRUE,
               gap = grid::unit(gap, "cm"),
               border = FALSE,
               gp = gp,
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







#' GSEA plot modified from 'enrichplot::gseaplot2'
#'
#' @param gsearesults
#' @param n
#'
#' @return
#' @export
#'
#' @examples
gseaplot <- function(gsearesults, n = 3){

  stopifnot(requireNamespace("enrichplot"))

  gsea_plot_modified <- function (x, geneSetID, title = "", color = "green", base_size = 20,
                                  rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE,
                                  ES_geom = "line"){

    ES_geom <- match.arg(ES_geom, c("line", "dot"))
    geneList <- position <- NULL
    if (length(geneSetID) == 1) {
      gsdata <- gsInfo(x, geneSetID)
    }
    else {
      gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
    }
    p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) +
      theme(panel.grid.major = element_blank(),
            axis.text = element_text(colour = "black"),
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
      scale_x_continuous(expand = c(0, 0))
    if (ES_geom == "line") {
      es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description),
                            size = 1)
    }
    else {
      es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description),
                             size = 1, data = subset(gsdata, position == 1))
    }
    p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8),
                                  text = element_text(size = 0.7 * base_size),
                                  legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
    p.res <- p.res + ylab("Enrichment Score") + theme(axis.text.x = element_blank(),
                                                      axis.ticks.x = element_blank(), axis.line.x = element_blank(),
                                                      plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2,
                                                                           unit = "cm"))
    i <- 0
    for (term in unique(gsdata$Description)) {
      idx <- which(gsdata$ymin != 0 & gsdata$Description ==
                     term)
      gsdata[idx, "ymin"] <- i
      gsdata[idx, "ymax"] <- i + 1
      i <- i + 1
    }
    p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin,
                                                             ymax = ~ymax, color = ~Description)) + xlab(NULL) +
      ylab(NULL) + theme_classic(base_size) + theme(legend.position = "none",
                                                    plot.margin = margin(t = -0.1, b = 0, unit = "cm"),
                                                    axis.ticks = element_blank(), axis.text = element_blank(),
                                                    axis.line.x = element_blank()) + scale_x_continuous(expand = c(0,
                                                                                                                   0)) + scale_y_continuous(expand = c(0, 0))
    if (length(geneSetID) == 1) {
      v <- seq(1, sum(gsdata$position), length.out = 9)
      inv <- findInterval(rev(cumsum(gsdata$position)), v)
      if (min(inv) == 0)
        inv <- inv + 1
      col = c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
      ymin <- min(p2$data$ymin)
      yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
      xmin <- which(!duplicated(inv))
      xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
      d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin,
                      xmax = xmax, col = col[unique(inv)])
      p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax,
                                ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d,
                           alpha = 0.9, inherit.aes = FALSE)
    }
    df2 <- p$data
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x,
                                               y = ~y, yend = 0), color = "grey")
    p.pos <- p.pos + ylab("Ranked list metric") + xlab("Rank in Ordered Dataset") +
      theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2,
                                 l = 0.2, unit = "cm"))
    if (!is.null(title) && !is.na(title) && title != "")
      p.res <- p.res + ggtitle(title)
    if (length(color) == length(geneSetID)) {
      p.res <- p.res + scale_color_manual(values = color) + theme(legend.position = "left", legend.direction = "vertical")
      if (length(color) == 1) {
        p.res <- p.res + theme(legend.position = "left", legend.direction = "vertical")
        p2 <- p2 + scale_color_manual(values = "black")
      }
      else {
        p2 <- p2 + scale_color_manual(values = color)
      }
    }
    if (pvalue_table) {
      pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
      pd <- pd[order(pd[, 1], decreasing = FALSE), ]
      rownames(pd) <- pd$Description
      pd <- pd[, -1]
      pd <- round(pd, 4)
      tp <- tableGrob2(pd, p.res)
      p.res <- p.res
    }

    p2 <- p2 + theme(text = element_text(size = base_size*0.7))

    p3 <- cowplot::get_legend(p.res)
    p.res <- p.res + theme(legend.position = "none")

    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(),
                                           axis.ticks.x = element_line(), axis.text.x = element_text())
    if (length(subplots) == 1)
      return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2,
                                                        r = 0.2, b = 0.2, l = 0.2, unit = "cm")))

    p <- plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
    plot_grid(p, p3, ncol = 1, rel_heights = c(1, 0.025 + 0.025*length(geneSetID)))
  }


  gg <- gsea_plot_modified(gsearesults, geneSetID = 1:n, base_size = 16, subplots = 1:2, rel_heights = c(1, 0.15),
                           color =  c("#eb0000", "#4d58ff", "#fac000", "#12c939", "#36a5c9", topo.colors(10))[1:n])
  gg <- gg +
    theme_bw() +
    theme(plot.title = element_text(size = 18, hjust=0.5),
          plot.margin = unit(c(1,2.5,1,2),"cm"),
          panel.grid = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())

  gg
}












#' GSEA bar plot
#'
#' @param results
#' @param x
#' @param y
#' @param label
#' @param sort_by
#' @param sort_abs
#' @param top_by
#' @param top_n
#' @param top_n_up
#' @param top_n_down
#' @param labsize
#' @param reverse
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ggseabar <- function(results, x = NULL, y = NULL, label = NULL, sort_by = NULL, sort_abs = FALSE, top_by = NULL, top_n = 20, top_n_up = NULL, top_n_down = NULL, labsize = 12, reverse = FALSE, ...){

  # Input
  results <- as.data.frame(results)

  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  label <- rlang::enquo(label)
  sort_by <- rlang::enquo(sort_by)
  top_by <- rlang::enquo(top_by)

  ordered_guess <- function(x, ref){
    guess <- grep(paste0(x, collapse = "|"), ref, value = TRUE)
    guess <- sort(factor(guess, ordered = TRUE, levels = x))
    guess[1]
  }

  if (rlang::quo_is_null(x)){
    x <- ordered_guess(c("padj","p.adjust","padjust","fdr", "FDR", "pvalue","pval"), colnames(results))
    x <- rlang::sym(as.character(x))
  }

  if (rlang::quo_is_null(y)){
    y <- ordered_guess( c("pathway", "Pathway", "description", "Description", "id", "ID", "geneset"), colnames(results))
    y <- rlang::sym(as.character(y))
  }

  if (rlang::quo_is_null(label)) label <- y
  if (rlang::quo_is_null(sort_by)) sort_by <- x
  if (rlang::quo_is_null(top_by)) top_by <- sort_by

  # Processing
  if (is.null(top_n_up)){
    results <- dplyr::arrange(results, dplyr::desc(abs(!!top_by)))
    results <- results[1:top_n,, drop = FALSE]
  } else {
    results <- dplyr::arrange(results, dplyr::desc(!!top_by))
    if (top_n_up > 0) ix.up <- 1:top_n_up else ix.up <- c()
    if (top_n_down > 0) ix.down <- (nrow(results) - top_n_down + 1):nrow(results) else ix.down <- c()
    ix <- unique(c(ix.up, ix.down))
    results <- results[ix,, drop = FALSE]
  }

  if (sort_abs == TRUE){
    results <- dplyr::arrange(results, dplyr::desc(abs(!!sort_by)))
  } else {
    results <- dplyr::arrange(results, dplyr::desc(!!sort_by))
  }

  if (reverse == TRUE){
    results[[rlang::as_name(y)]] <- factor(results[[rlang::as_name(y)]], ordered = TRUE, levels = results[[rlang::as_name(y)]][!duplicated(results[[rlang::as_name(y)]])])
  } else {
    results[[rlang::as_name(y)]] <- factor(results[[rlang::as_name(y)]], ordered = TRUE, levels = rev(results[[rlang::as_name(y)]][!duplicated(results[[rlang::as_name(y)]])]))
  }


  gg <- ggplot2::ggplot(results, aes(x = !!x, y = !!y, label = !!label, ...)) + ggplot2::geom_bar(stat = "identity")
  gg <- gg + ggplot2::geom_vline(xintercept = 0, color = "grey40")

  res_left <- subset(results, results[[rlang::as_name(x)]] < 0)
  if (nrow(res_left) > 0) gg <- gg + ggplot2::geom_text(size = labsize/ggplot2::.pt, data = res_left, nudge_x = 0.25, hjust = 0)
  res_right <- subset(results, results[[rlang::as_name(x)]] > 0)
  if (nrow(res_right) > 0) gg <- gg + ggplot2::geom_text(size = labsize/ggplot2::.pt, data = res_right, nudge_x = -0.25, hjust = 1)

  lims <- ceiling(abs(range(results[[rlang::as_name(x)]])))
  gg <- gg + ggplot2::ylab("")
  gg <- gg + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                            panel.border = ggplot2::element_blank(),
                            plot.background = ggplot2::element_blank(),
                            axis.ticks.y = ggplot2::element_blank(),
                            plot.title = ggplot2::element_text(hjust = 0.5))


  gg <- gg + ggplot2::xlim(c(-max(abs(lims)), max(abs(lims))))
  gg <- gg + ggplot2::scale_fill_gradient2(low = "#0e7aed", mid = "white", high = "#db3b25", midpoint = 0, oob = scales::squish)

  gg
}
















