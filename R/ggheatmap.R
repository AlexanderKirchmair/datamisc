


# add significance labels
# font sizes
# zscore
# add column names as separate plots? (should align with dendro)


# # Example data:
# data <- rmat(80, 30)
# data <- data-mean(data)
# rowdata <- data.frame(row.names = rownames(data), group = sample(c("groupA", "groupB"), nrow(data), replace = TRUE))
# rowdata$group2 <- paste0("G", sample(1:3, nrow(data), replace = TRUE))
# coldata <- data.frame(row.names = colnames(data), class = sample(c("A", "B", "C"), ncol(data), replace = TRUE))
#
#
# ggheatmap(data)
#
# ggheatmap(data, rowdata = rowdata, rownames = T, rowdata_names = T, check = F,widths = c(1,2,1),coldata = coldata,
#           layout = list("top" = "col_dendro", "bottom" = c("col_anno"), "left" = c("row_dendro", "row_anno", "legend")), vjust = -1.5)
#
# hm <- ggheatmap(data, lab = "score", rowdata = rowdata, rownames = T, rowdata_names = F, check = F, plot = F)
# hm





ggheatmap <- function(data,
                      lab = "value", legend = TRUE,
                      rownames = TRUE, rowdata = NULL, rowdata_names = FALSE,
                      colnames = TRUE, coldata = NULL,
                      colors = NULL,
                      low = "#004cff", mid = "white", high = "#ffae00",
                      layout = NULL, widths = c(0.1,1,0.1), heights = c(0.1,1,0.1),
                      tidy = FALSE, check = FALSE, plot = TRUE,
                      cluster_rows = NULL, cluster_cols = NULL, dendro_rows = TRUE, dendro_cols = TRUE, cluster_NA = NULL, clusterFUN = NULL,
                      ...){




  ### Input arguments

  # determine whether to perform clustering based on number of rows/columns
  if (is.null(cluster_rows) & nrow(data) < 1000) cluster_rows <- TRUE
  if (is.null(cluster_cols) & ncol(data) < 1000) cluster_cols <- TRUE


  if (is.null(layout)){ layout <- list("top" = c("col_dendro"), "right" = c("legend"), "bottom" = c("col_anno"), "left" = c("row_dendro", "row_anno")) }


  gg <- list()
  row_clusters <- NULL
  col_clusters <- NULL





  ### Data formatting ----
  # both data and data_tidy are needed

  if (tidy == FALSE){

    # if input data are not tidy, convert to long format (used e.g. for ggplot2):
    data_tidy <- as.data.frame.table(data)
    colnames(data_tidy) <- c("y", "x", "value")
    data_tidy$x <- factor(data_tidy$x, ordered = TRUE, levels = rev(colnames(data))) # set cluster order
    data_tidy$y <- factor(data_tidy$y, ordered = TRUE, levels = rev(rownames(data))) # set cluster order

  } else {

    # if input data are tidy, convert to matrix format (used e.g. for clustering)
    data_tidy <- data
    data <- tidyr::spread(data_tidy, key = "x", value = "value")
    data <- tibble::column_to_rownames(data, var = "y")

  }


  ### Clustering ----
  # todo: update dendsort

  # substitution of NA values for clustering
  data_clust <- data
  if (!is.null(cluster_NA)){
    data_clust[is.na(data_clust)] <- cluster_NA
  }

  # clustering function (can be user-supplied)
  if (is.null(clusterFUN)){
    clusterFUN <- function(data){
      dendsort::dendsort(stats::hclust(stats::dist(data)))
    }
  }

  # Row clusters
  if (cluster_rows == TRUE){
    row_clusters <- clusterFUN(data)
    data_tidy$y <- factor(data_tidy$y, ordered = TRUE, levels = row_clusters$labels[row_clusters$order]) # reorder data according to the clustering
    if (dendro_rows == TRUE) gg$row_dendro <- .getDendrogram(cluster = row_clusters, ddata = ddata, side = "x", check = check)
  }

  # Column clusters
  if (cluster_cols == TRUE){
    col_clusters <- clusterFUN(t(data))
    data_tidy$x <- factor(data_tidy$x, ordered = TRUE, levels = col_clusters$labels[col_clusters$order]) # reorder data according to the clustering
    if (dendro_cols == TRUE) gg$col_dendro <- .getDendrogram(cluster = col_clusters, ddata = ddata, side = "y", check = check)
  }



  ### Colors ----

  if (is.null(colors)) colors <- list()

  colors <- c(colors,
              getColors(coldata[,!colnames(coldata) %in% names(colors), drop = FALSE]),
              getColors(rowdata[,!colnames(rowdata) %in% names(colors), drop = FALSE]))

  if (is.null(colors$hm)){
    if (any(data_tidy$value < 0)){
      colors$hm <- ggplot2::scale_fill_gradient2(name = lab, low = low, mid = mid, high = high)
    } else {
      colors$hm <- ggplot2::scale_fill_gradient(name = lab, low = low, high = high)
    }
  }



  ### Annotations ----

  ggrow <- .ggAnno(rowdata, type = "row", anno_names = rowdata_names, layout = layout, colors = colors, gg = gg, check = check, tidy = tidy, data_tidy = data_tidy, ...)
  gg <- c(gg, ggrow)

  ggcol <- .ggAnno(coldata, type = "col", anno_names = coldata_names, layout = layout, colors = colors, gg = gg, check = check, tidy = tidy, data_tidy = data_tidy, ...)
  gg <- c(gg, ggcol)




  ### Main heatmap ----

  gg$hm <- ggplot2::ggplot(data = data_tidy, mapping = ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.title = ggplot2::element_blank()) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::scale_y_discrete(position = "right") +
    colors$hm


  if (rownames == TRUE) gg$hm <- gg$hm + ggplot2::theme(axis.text.y = ggplot2::element_text())
  if (colnames == TRUE) gg$hm <- gg$hm + ggplot2::theme(axis.text.x = ggplot2::element_text())





  ### Combine plots ----

  # set shared theme elements
  gg <- lapply(gg, function(tmp) tmp + ggplot2::coord_cartesian(clip = "off"))
  gg <- lapply(gg, function(tmp) tmp + ggplot2::theme(axis.title = ggplot2::element_blank(),
                                                      plot.margin = ggplot2::margin(0,0,0,0),
                                                      panel.spacing = ggplot2::margin(0,0,0,0)) )

  if (check == TRUE){
    gg <- lapply(gg, function(tmp) tmp + ggplot2::theme(plot.background = ggplot2::element_rect(fill = rgb(0.9,0.9,0.9)),
                                                        panel.background = ggplot2::element_rect(fill = rgb(0.8,0.8,1)) ))
  }


  if ("legend" %in% unlist(layout)) gg$legend <- patchwork::guide_area()

  lix <- grepl("anno", names(gg))
  if (legend == FALSE) gg[lix] <- lapply(gg[lix], function(tmp) tmp + ggplot2::theme(legend.position = "none") )


  # arrange plots
  p <- .getLayout(gg, layout, widths, heights)
  p$cluster$row <- row_clusters
  p$cluster$col <- col_clusters



  ### Output ----

  class(p) <- c("ggheatmap", "gg", "list")
  if (plot == TRUE) p <- plot(p)
  p


}







greps <- function(patterns, x, ..., FUN = grepl){
  res <- sapply(patterns, function(p) FUN(x = x, pattern = p) )
  apply(res, 1, function(tmp) Reduce("|", tmp) )
}




.plotGgheatmap <- function(p, ...){

  patchwork::wrap_plots(p$plotlist,
                        guides = "collect",
                        design = p$layout,
                        widths = p$widths_adj,
                        heights = p$heights_adj,
                        ...)

}


plot.ggheatmap <- function(x, ...){
  .plotGgheatmap(p = x, ...)
}


print.ggheatmap <- function(x, ...){
  print(plot(x, ...))
}


.getLayout <- function(gg, layout, widths, heights){

  # arrange the plotlist into top, right, bottom, left parts based on layout
  g <- lapply(layout, function(ltmp) gg[greps(ltmp, names(gg))] )

  if ("legend" %in% names(g$left)){
    ix <- grep("legend", names(g$left))
    g$left <- c(g$left[ix], g$left[-ix])
  }

  # add main
  g$main <- list("A" = gg$hm)
  if (!is.null(g$top)) names(g$top) <- setdiff(LETTERS, "A")[seq(g$top)]
  if (!is.null(g$left)) names(g$left) <- setdiff(LETTERS, c("A", names(g$top)))[seq(g$left)]
  if (!is.null(g$right)) names(g$right) <- setdiff(LETTERS, c("A", names(g$top), names(g$left)))[seq(g$right)]
  if (!is.null(g$bottom)) names(g$bottom) <- setdiff(LETTERS, c("A", names(g$top), names(g$left), names(g$right)))[seq(g$bottom)]

  col <- c(names(g$top), names(g$main), names(g$bottom))
  row <- c(names(g$left), names(g$main), names(g$right))
  mat <- matrix(nrow = length(col), ncol = length(row))

  # add corners (annotation titles)

  mat[which(col == "A"),] <- row
  mat[,which(row == "A")] <- col
  mat[is.na(mat)] <- "#"


  mat <- paste0(apply(mat, 1, function(x) paste0(c(x, "\n"), collapse = "") ), collapse = " ")

  p <- unlist(setNames(g, NULL), recursive = FALSE)


  widths_adj <- c(rep(widths[1]/length(g$left), length(g$left)),
                  widths[2],
                  rep(widths[3]/length(g$right), length(g$right)))

  heights_adj <- c(rep(heights[1]/length(g$top), length(g$top)),
                   heights[2],
                   rep(heights[3]/length(g$bottom), length(g$bottom)))


  list(plotlist = p, layout = mat, widths_adj = widths_adj, heights_adj = heights_adj)
}





.ggAnno <- function(annodata, type, layout, colors, gg, data_tidy, anno_names = FALSE, check = FALSE, tidy = FALSE, ...){

  if (is.null(annodata)) return(NULL)
  stopifnot(type %in% c("row", "col"))


  if ("gg" %in% class(annodata)){
    return(annodata)

  } else {

    annodata_tidy <- NULL

    if (tidy == FALSE){

      annodata_tidy <- as.data.frame(tidyr::pivot_longer(dplyr::mutate(annodata, id = rownames(annodata)), -id))
      annodata_tidy$name <- factor(annodata_tidy$name, ordered = TRUE, levels = rev(colnames(annodata))) # set cluster order

      if (type == "row"){
        #colnames(annodata_tidy) <- c("y", "x", "value")
        annodata_tidy$id <- factor(annodata_tidy$id, ordered = TRUE, levels = levels(data_tidy$y)) # set cluster order
      }

      if (type == "col"){
        #colnames(annodata_tidy) <- c("x", "y", "value")
        annodata_tidy$id <- factor(annodata_tidy$id, ordered = TRUE, levels = levels(data_tidy$x)) # set cluster order
      }


    } else {

      # pull from tidy data








    }

    if (!is.null(annodata_tidy)){

      annodata_list <- setNames(split(annodata_tidy, annodata_tidy$name), NULL)
      names(annodata_list) <- sapply(annodata_list, function(x) unique(x$name))

      if (type == "row"){
        annodata_list <- lapply(annodata_list, function(tmp) dplyr::rename(.data = tmp, x = name, y = id))
      } else if (type == "col"){
        annodata_list <- lapply(annodata_list, function(tmp) dplyr::rename(.data = tmp, x = id, y = name))
      }



      # annotation plot (should have the same y axis as the main heatmap)
      anno_plots <- lapply(setNames(names(annodata_list), names(annodata_list)), function(tmp){

        datatmp <- annodata_list[[tmp]]

        gg <- ggplot2::ggplot(data = datatmp, mapping = ggplot2::aes(x = x, y = y, fill = value)) +
          ggplot2::theme_void() +
          ggplot2::geom_tile() +
          ggplot2::scale_x_discrete(expand = c(0,0))

        if (type == "row" & any(grepl("row_anno", layout$right))) gg <- gg + ggplot2::scale_y_discrete(position = "right")
        if (type == "col" & any(grepl("col_anno", layout$top))) gg <- gg + ggplot2::scale_x_discrete(position = "top")

        # add labels outside of the plotting area (to avoid gaps between heatmaps on top/bottom annotations)
        if (type == "row"){
          gg <- gg + ggplot2::annotate("text",
                                       x = tmp,
                                       y = rev(levels(datatmp$y))[1],
                                       label = tmp,
                                       vjust = -1.5)
        } else {
          gg <- gg + ggplot2::annotate("text",
                                       x = rev(levels(datatmp$x))[1],
                                       y = tmp,
                                       label = tmp,
                                       hjust = -1)
        }



        # anno_titles
        # ggplot2::ggplot(data.frame(x = tmp, y = 1), aes(x, y, label = x)) + theme_void() + geom_text()


        if ("Scale" %in% class(colors[[tmp]])){
          gg <- gg + colors[[tmp]]

        } else {

          ## add scale for numerical data #########



          gg <- gg + ggplot2::scale_fill_manual(name = tmp, values = colors[[tmp]])
        }

        gg
      })




      if (type == "row"){
        if (check == TRUE | (anno_names == TRUE & any(grepl("anno", layout$right)))) anno_plots[[length(anno_plots)]] <- anno_plots[[length(anno_plots)]] + ggplot2::theme(axis.text.y = ggplot2::element_text())
        if (check == TRUE | (anno_names == TRUE & any(grepl("anno", layout$left)))) anno_plots[[1]] <- anno_plots[[1]] + ggplot2::theme(axis.text.y = ggplot2::element_text())
      } else if (type == "col"){


      }
    }

    if (type == "row") anno_plots <- setNames(anno_plots, paste0("row_anno_", names(anno_plots)))
    if (type == "col") anno_plots <- setNames(anno_plots, paste0("col_anno_", names(anno_plots)))

    anno_plots

  }

}








.getDendrogram <- function(cluster, ddata, side = "x", check = FALSE){

  ddata <- ggdendro::dendro_data(as.dendrogram(cluster), type = "rectangle")
  ddata$labels$label <- factor(ddata$labels$label, ordered = TRUE, levels = ddata$labels$label[ddata$labels$x])

  dendro <- ggplot2::ggplot() + cowplot::theme_nothing()

  if (side == "x"){
    dendro <- dendro + ggplot2::geom_text(data = ddata$labels, color = ifelse(check == TRUE, "blue", NA), hjust = 1, vjust = 0.5, mapping = ggplot2::aes(x = 0, y = label, label = label)) + # this adds the correct x scale
      ggplot2::geom_segment(data = ddata$segments, ggplot2::aes(x = -y, y = x, xend = -yend, yend = xend)) +
      ggplot2::scale_x_continuous(expand = c(0,0))
  } else {
    dendro <- dendro + ggplot2::geom_text(data = ddata$labels, color = ifelse(check == TRUE, "blue", NA),  hjust = 0.5, vjust = 0, mapping = ggplot2::aes(x = label, y = 0, label = label)) + # this adds the correct y scale
      ggplot2::geom_segment(data = ddata$segments, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
      ggplot2::scale_y_continuous(expand = c(0,0))

  }

  dendro
}


