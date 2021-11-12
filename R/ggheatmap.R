


# add significance labels



# # Example data:
# data <- rmat(100, 30)
# rowdata <- data.frame(row.names = rownames(data), group = sample(c("groupA", "groupB"), nrow(data), replace = TRUE))
# rowdata$group2 <- paste0("G", sample(1:3, nrow(data), replace = TRUE))
# coldata <- data.frame(row.names = colnames(data), class = sample(c("A", "B", "C"), ncol(data), replace = TRUE))
#
#
# ggheatmap(data-mean(data), rowdata = rowdata, rownames = F, rowdata_names = T, check = F,widths = c(1,2,1),
#           layout = list("top" = "col_dendro", "bottom" = c("col_anno"), "left" = c("row_dendro", "row_anno")))
#
#
# ggheatmap(data-mean(data), rowdata = rowdata, rownames = F, rowdata_names = T, check = T)
#
# gg <- ggheatmap(data-mean(data), rowdata = rowdata, rownames = F, rowdata_names = T, check = F, return_list = T)
#
# # scale_fill for annotations: either 1 or all cols or combining multiple individual ggplots
#





ggheatmap <- function(data,
                      rownames = TRUE, rowdata = NULL, rowdata_names = FALSE,
                      colnames = TRUE, coldata = NULL,
                      colors = NULL,
                      layout = NULL, widths = c(0.1,1,0.1), heights = c(0.1,1,0.1),
                      tidy = FALSE, check = FALSE,
                      cluster_rows = NULL, cluster_cols = NULL, dendro_rows = TRUE, dendro_cols = TRUE, cluster_NA = NULL, clusterFUN = NULL,
                      return_list = FALSE, ...){



  gg <- list()
  if (is.null(colors)) colors <- list()


  if (is.null(layout)){ layout <- list("top" = c("col_dendro"), "right" = c("row_anno"), "bottom" = c("col_anno"), "left" = c("row_dendro")) }




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


    rowdata_list <- setNames(split(rowdata_tidy, rowdata_tidy$x), NULL)
    names(rowdata_list) <- sapply(rowdata_list, function(x) unique(x$x))


    # annotation plot (should have the same y axis as the main heatmap)
    row_anno_plots <- lapply(setNames(names(rowdata_list), names(rowdata_list)), function(tmp){

      rowtmp <- rowdata_list[[]]

      gg <- ggplot2::ggplot(data = rowtmp, mapping = ggplot2::aes(x = x, y = y, fill = value)) +
        ggplot2::theme_void() +
        ggplot2::geom_tile() +
        ggplot2::scale_x_discrete(expand = c(0,0)) +
        ggplot2::scale_y_discrete(position = "right")

      # add labels outside of the plotting area (to avoid gaps between heatmaps on top/bottom annotations)
      x <- setdiff(unique(rowtmp$x), "id")
      gg + ggplot2::annotate("text",
                              x = x,
                              y = rev(levels(rowtmp$y))[1],
                              label = x,
                              vjust = -1.5)

    })


    # gg$row_anno + ggplot2::facet_wrap(~ x, scales = "free_x")
    # $row_anno + ggplot2::guides(fill = ggplot2::guide_legend(ncol= ncol(rowdata)))

    if (check == TRUE | rowdata_names == TRUE) row_anno_plots[[length(row_anno_plots)]] <- row_anno_plots[[length(row_anno_plots)]] + ggplot2::theme(axis.text.y = ggplot2::element_text())

    gg <- c(gg, setNames(row_anno_plots, paste0("row_anno_", names(row_anno_plots))))



  }













  ### Main heatmap ----

  gg$hm <- ggplot2::ggplot(data = data_tidy, mapping = ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::theme_void() +
    ggplot2::theme(axis.title = ggplot2::element_blank()) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_x_discrete(position = "bottom") +
    ggplot2::scale_y_discrete(position = "right")


  if (any(data_tidy$value < 0)){
    gg$hm <-  gg$hm + ggplot2::scale_fill_gradient2(low = "#004cff", mid = "white", high = "#ffae00")
  } else {
    gg$hm <-  gg$hm + ggplot2::scale_fill_gradient(low = "#ededed", high = "#004cff")
  }



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

  p <- .getLayout(gg, layout, widths, heights)

  gg_comb <- patchwork::wrap_plots(p$plotlist,
                                   guides = "collect",
                                   design = p$layout,
                                   widths = p$widths_adj,
                                   heights = p$heights_adj)

  # row/col-titles
  # overall title

  gg_comb

}


greps <- function(patterns, x, ..., FUN = grepl){
  res <- sapply(patterns, function(p) FUN(x = x, pattern = p) )
  apply(res, 1, function(tmp) Reduce("|", tmp) )
}




.getLayout <- function(gg, layout, widths, heights){

  g <- lapply(layout, function(ltmp) gg[greps(ltmp, names(gg))] )

  g$main <- list("A" = gg$hm)
  if (!is.null(g$top)) names(g$top) <- setdiff(LETTERS, "A")[seq(g$top)]
  if (!is.null(g$left)) names(g$left) <- setdiff(LETTERS, c("A", names(g$top)))[seq(g$left)]
  if (!is.null(g$right)) names(g$right) <- setdiff(LETTERS, c("A", names(g$top), names(g$left)))[seq(g$right)]
  if (!is.null(g$bottom)) names(g$bottom) <- setdiff(LETTERS, c("A", names(g$top), names(g$left), names(g$right)))[seq(g$bottom)]

  col <- c(names(g$top), names(g$main), names(g$bottom))
  row <- c(names(g$left), names(g$main), names(g$right))
  mat <- matrix(nrow = length(col), ncol = length(row))

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

  # mat <- matrix(0, nrow = 3, ncol = 3)
  # mat[2,2] <- NA
  # mat[2,1] <- length(g$left)
  # mat[1,2] <- length(g$top)
  # mat[2,3] <- length(g$right)
  # mat[3,2] <- length(g$bottom)

#
#   # generate layout - columns
#   i <- 1
#   n <- ncol(mat)
#   while(i <= n){
#     n <- ncol(mat)
#     cent <- which(is.na(rowSums(mat)))
#     if (sum(mat[cent,i], na.rm = TRUE) > 1){
#       if (i == 1) tmp0 <- NULL else tmp0 <- mat[,1:(i-1)]
#       tmp <- matrix(rep(mat[,i], mat[cent,i]), nrow = nrow(mat))
#       tmp[cent,] <- 1
#       if (i < ncol(mat)-1) tmp2 <- mat[,(i+1):ncol(mat)] else tmp2 <- NULL
#       mat <- cbind(tmp0, tmp, tmp2)
#     }
#     i <- i+1
#   }
#
#
#   # generate layout - rows
#   i <- 1
#   n <- nrow(mat)
#   while(i <= n){
#     n <- nrow(mat)
#     cent <- which(is.na(colSums(mat)))
#     if (sum(mat[i,cent], na.rm = TRUE) > 1){
#       if (i == 1) tmp0 <- NULL else tmp0 <- mat[1:(i-1),]
#       tmp <- matrix(rep(mat[i,], mat[i,cent]), nrow = nrow(mat), byrow = TRUE)
#       tmp[,cent] <- 1
#       if (i < nrow(mat)-1) tmp2 <- mat[(i+1):nrow(mat),] else tmp2 <- NULL
#       mat <- rbind(tmp0, tmp, tmp2)
#     }
#     i <- i+1
#   }
#
#
#   mat <- mat[rowSums(mat) != 0 | is.na(rowSums(mat)), colSums(mat) != 0 | is.na(colSums(mat))]
#
#   lmat <- mat * NA
#   lmat[is.na(mat)] <- "Z"
#   lmat[mat == 0] <- "#"
#   for (i in 1:nrow(mat)){
#     for (j in 1:ncol(mat)){
#       if (is.na(lmat[i,j])) lmat[i,j] <- setdiff(LETTERS, unique(as.vector(lmat)))[1]
#     }
#   }
#
#   # make named plotlist
#   ind <- which(lmat == "Z", arr.ind = TRUE)
#
#
#   top <- lmat[,ind[,"col"]]


  list(plotlist = p, layout = mat, widths_adj = widths_adj, heights_adj = heights_adj)
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

