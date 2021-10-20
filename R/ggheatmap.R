




# # Example data:
# data <- rmat(100, 30)
# rowdata <- data.frame(row.names = rownames(data), group = sample(c("groupA", "groupB"), nrow(data), replace = TRUE))
# rowdata$group2 <- paste0("G", sample(1:3, nrow(data), replace = TRUE))
# coldata <- data.frame(row.names = colnames(data), class = sample(c("A", "B", "C"), ncol(data), replace = TRUE))
#
#
# ggheatmap(data, rowdata = rowdata, rownames = F, rowdata_names = T, check = F)




ggheatmap <- function(data,
                      rownames = TRUE, rowdata = NULL, rowdata_names = FALSE,
                      colnames = TRUE, coldata = NULL,
                      tidy = FALSE, check = FALSE,
                      cluster_rows = NULL, cluster_cols = NULL, dendro_rows = TRUE, dendro_cols = TRUE, cluster_NA = NULL, clusterFUN = NULL,
                      return_list = FALSE, ...){



  gg <- list()


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


  # determine whether to perform clustering based on number of rows/columns
  if (is.null(cluster_rows)) cluster_rows <- TRUE
  if (is.null(cluster_cols)) cluster_cols <- TRUE

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

    if (dendro_rows == TRUE){ # plot if TRUE

      ddata <- ggdendro::dendro_data(as.dendrogram(row_clusters), type = "rectangle")
      ddata$labels$label <- factor(ddata$labels$label, ordered = TRUE, levels = ddata$labels$label[ddata$labels$x])

      gg$row_dendro <- ggplot2::ggplot() +
        cowplot::theme_nothing() +
        ggplot2::geom_text(data = ddata$labels, color = ifelse(check == TRUE, "blue", NA), hjust = 1, vjust = 0.5, mapping = ggplot2::aes(x = 0, y = label, label = label)) + # this adds the correct x scale
        ggplot2::geom_segment(data = ddata$segments, ggplot2::aes(x = -y, y = x, xend = -yend, yend = xend)) +
        ggplot2::scale_x_continuous(expand = c(0,0))

    }
  }

  # Column clusters
  if (cluster_cols == TRUE){

    col_clusters <- clusterFUN(t(data))
    data_tidy$x <- factor(data_tidy$x, ordered = TRUE, levels = col_clusters$labels[col_clusters$order]) # reorder data according to the clustering

    if (dendro_cols == TRUE){ # plot if TRUE

      ddata <- ggdendro::dendro_data(as.dendrogram(col_clusters), type = "rectangle")
      ddata$labels$label <- factor(ddata$labels$label, ordered = TRUE, levels = ddata$labels$label[ddata$labels$x])

      gg$col_dendro <- ggplot2::ggplot() +
        cowplot::theme_nothing() +
        ggplot2::geom_text(data = ddata$labels, color = ifelse(check == TRUE, "blue", NA),  hjust = 0.5, vjust = 0, mapping = ggplot2::aes(x = label, y = 0, label = label)) + # this adds the correct y scale
        ggplot2::geom_segment(data = ddata$segments, ggplot2::aes(x = x, y = y, xend = xend, yend = yend)) +
        ggplot2::scale_y_continuous(expand = c(0,0))

    }
  }




  ### Annotations ----




  rowdata_tidy <- NULL

  if (tidy == FALSE){

    rowdata_tidy <- as.data.frame(tidyr::pivot_longer(dplyr::mutate(rowdata, id = rownames(rowdata)), -id))
    colnames(rowdata_tidy) <- c("y", "x", "value")
    rowdata_tidy$x <- factor(rowdata_tidy$x, ordered = TRUE, levels = rev(colnames(rowdata))) # set cluster order
    rowdata_tidy$y <- factor(rowdata_tidy$y, ordered = TRUE, levels = levels(data_tidy$y)) # set cluster order

  } else {

    # pull from tidy data

  }






  # Row annotation

  if (!is.null(rowdata_tidy)){

    # annotation plot (should have the same y axis as the main heatmap)
    gg$row_anno <- ggplot2::ggplot(data = rowdata_tidy, mapping = ggplot2::aes(x = x, y = y, fill = value)) +
      ggplot2::theme_void() +
      ggplot2::geom_tile() +
      ggplot2::scale_x_discrete(expand = c(0,0)) +
      ggplot2::scale_y_discrete(position = "right")

    # add labels outside of the plotting area (to avoid gaps between heatmaps on top/bottom annotations)
    x <- setdiff(levels(rowdata_tidy$x), "id")
    gg$row_anno <- gg$row_anno + ggplot2::annotate("text",
                                    x = x,
                                    y = rev(levels(rowdata_tidy$y))[1],
                                    label = x,
                                    vjust = -1.5)

    if (check == TRUE | rowdata_names == TRUE) gg$row_anno <- gg$row_anno + ggplot2::theme(axis.text.y = ggplot2::element_text())

  }













  ### Main heatmap ----

  gg$hm <- ggplot2::ggplot(data = data_tidy, mapping = ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.title = ggplot2::element_blank()) +
    ggplot2::geom_tile() +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::scale_y_discrete(position = "right")

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


  if (return_list == TRUE){ # return list of plots instead of the combined final plot
    # combine legends from all plots
    # legends <- lapply(gg, cowplot::get_legend)
    # legends <- legends[lengths(legends) > 0]
    # gg <- lapply(gg, function(tmp) tmp + ggplot2::theme(legend.position = "none"))
    # gg$legends <- patchwork::wrap_plots(legends, ncol = ceiling(length(legends)/2))
    return(gg)
  }

  # arrange plots

  p <- list()
  p$A <- gg$col_dendro
  p$B <- gg$row_dendro
  p$C <- gg$hm
  p$D <- gg$row_anno
  p$E <- NULL


  layout <- "#A#
             BCD
             #E#"

  for (i in setdiff(LETTERS, names(p))){
    layout <- gsub(i, "#", layout)
  }

  gg_comb <- patchwork::wrap_plots(p,
                                   guides = "collect",
                                   design = layout,
                                   widths = c(1,3,1),
                                   heights = c(1,3,1))
  gg_comb

}


# # Customize:
#
#
# gh <- ggheatmap(data, rowdata = rowdata, rownames = F, rowdata_names = T, check = F)
# gh[[3]] <- gh[[3]] + ggplot2::scale_fill_gradient2(low = "red", mid = "white", high = "blue")
# gh
#
#
# ggl <- ggheatmap(data, rowdata = rowdata, rownames = T, rowdata_names = T, check = F, return_list = TRUE)
#
# names(ggl)
# ggl2 <- setNames(ggl, c("d", "e", "r", "H"))
#
# layout <- "##e
#            drH
#            ###"
#
# patchwork::wrap_plots(ggl2,
#                       guides = "collect",
#                       design = layout,
#                       widths = c(1,1,3),
#                       heights = c(1,3,1))

