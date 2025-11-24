

#' Volcano plots using ggplot2
#'
#' @importFrom magrittr %<>%
#' @importFrom ggplot2 %+%
#' @param data Dataframe
#' @param x Variable to plot on x-axis
#' @param y Variable to plot on y-axis
#' @param color
#' @param label
#' @param shape
#' @param nlabels
#' @param lab_size
#' @param repel
#' @param attract
#' @param box.padding
#' @param max_overlaps
#' @param seed
#' @param ptres
#' @param clip
#' @param symlim
#' @param expand
#' @param nbreaks_x
#' @param nbreaks_y
#' @param color_up
#' @param color_down
#' @param color_nonsig
#' @param title
#' @param title_size
#' @param point_size
#' @param scale_size
#' @param axis_size
#' @param leg_size
#' @param lwd
#' @param at_zero
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' data.frame(row.names = LETTERS, lfc = runif(length(LETTERS), -3, 3), padj = runif(length(LETTERS))) |> ggvolcano()
ggvolcano <- function(data, x = NULL, y = NULL, color = NULL, label = NULL, shape = NULL, stroke = NA,
                      nlabels = NULL, lab_size = 12, labface = "plain", repel = 1.5, attract = NULL, box.padding = 0.5, max_overlaps = Inf, seed = 123,
                      ptres = 0.05, clip = FALSE, symlim = TRUE, expand = c(0,0), nbreaks_x = 7, nbreaks_y = 7,
                      xlim = NULL, ylim = NULL,
                      color_up = "#eb9d0e", color_down = "#146bc7", color_nonsig = "#4d4d4d",
                      label_up = "up", label_down = "down", label_nonsig = "not signif.", show_nonsig = TRUE, show_grid = TRUE,
                      title = NULL, title_size = NULL, point_size = 2, scale_size = FALSE, axis_size = NULL, leg_size = NULL, leg_key_size=3.5,
                      lwd = 0.8, at_zero = FALSE, clip_frame = "off", ...){

  ### Function to plot volcano plots using ggplot.


  # PARSE INPUTS

  data <- as.data.frame(data)

  x <- rlang::enquo(x)
  y <- rlang::enquo(y)

  if (rlang::quo_is_null(x)){ x <- rlang::sym(grep("lfc|log2FoldChange|logFC|log2FC|nes", names(data), value = TRUE, ignore.case = TRUE)[1]) }
  if (rlang::quo_is_null(y)){ y <- rlang::sym(grep("padj|fdr", names(data), value = TRUE, ignore.case = TRUE)[1]) }

  data$x <- data[[rlang::as_name(x)]]
  data$y <- -log10(data[[rlang::as_name(y)]])
  data <- data[!is.na(data$x) & !is.na(data$y),]

  data$xtmp <- data$x
  data$xtmp[is.infinite(data$xtmp)] <- max(abs(data$xtmp[!is.infinite(data$xtmp)])) * sign(data$xtmp[is.infinite(data$xtmp)])
  data$ytmp <- data$y
  data$ytmp[is.infinite(data$ytmp)] <- max(abs(data$ytmp[!is.infinite(data$ytmp)])) * sign(data$ytmp[is.infinite(data$ytmp)])
  data$score <- abs(as.numeric(scale(data$xtmp, center = FALSE))) + abs(as.numeric(scale(data$ytmp, center = FALSE)))
  data$score[is.na(data$score)] <- 0

  data$class <- label_nonsig
  data$class[data[[rlang::as_name(y)]] <= ptres & data$x > 0] <- label_up
  data$class[data[[rlang::as_name(y)]] <= ptres & data$x < 0] <- label_down

  data$score[data$class == label_nonsig] <- data$score[data$class == label_nonsig] * 0.001

  if (is.null(title_size)) title_size <- lab_size
  if (is.null(axis_size)) axis_size <- lab_size
  if (is.null(leg_size)) leg_size <- lab_size


  shape <- rlang::enquo(shape)


  # LABELS

  label <- rlang::enquo(label)

  if (rlang::quo_is_null(label)){
    data[["label"]] <- rownames(data)
  } else {
    data[["label"]] <- data[[rlang::as_name(label)]]
  }

  data <- data[order(data$score, decreasing = TRUE),]

  if (is.null(nlabels)){ nlabels <- min(20, ceiling(nrow(data)/10)) }
  if (is.infinite(nlabels)){ nlabels <- nrow(data) }
  data$do_label <- FALSE
  nlabels_left <- nlabels_right <- 0
  if (nrow(subset(data, x < 0)) > 0) nlabels_left <- ceiling(nlabels/2 * max(subset(data, x < 0)$score, na.rm = TRUE) / max(data$score, na.rm = TRUE))
  if (nrow(subset(data, x > 0)) > 0) nlabels_right <- ceiling(nlabels/2 * max(subset(data, x > 0)$score, na.rm = TRUE) / max(data$score, na.rm = TRUE))
  if (is.na(nlabels_left)) nlabels_left <- 0
  if (is.na(nlabels_right)) nlabels_right <- 0

  data$do_label[which(data$x < 0)[1:nlabels_left]] <- TRUE
  data$do_label[which(data$x > 0)[1:nlabels_right]] <- TRUE
  data$do_label[is.na(data$do_label)] <- FALSE

  if (sum(data$do_label) < nlabels){ data$do_label[!data$do_label][1:(nlabels-sum(data$do_label))] <- TRUE }
  data$label[!data$do_label] <- ""
  data$do_label[data$label == ""] <- FALSE

  data$do_label[data$class == label_nonsig] <- FALSE

  # COLORS

  color <- rlang::enquo(color)
  color_user_def <- rlang::quo_is_null(color)
  if (color_user_def){
    color <- rlang::sym("class")
  }


  # LIMITS

  sigdata <- subset(data, class != label_nonsig)
  xylimits <- list(xlim = getLimits(sigdata$xtmp, clip = clip, expand = expand[1]), ylim = getLimits(sigdata$ytmp, clip = clip, expand = expand[2], negative = FALSE))
  if (symlim == TRUE){ xylimits$xlim <- c("min" = -max(abs(xylimits$xlim)), "max" = max(abs(xylimits$xlim))) }
  data$xorg <- data$x
  data$yorg <- data$y

  if (!is.null(xlim)) xylimits$xlim <- setNames(xlim, c("min", "max"))
  if (!is.null(ylim)) xylimits$ylim <- setNames(ylim, c("min", "max"))

  # set limits and axis tick labels
  # if clip==TRUE and any points are cut off, and recalculate limits based on breaks
  # if clip==TRUE and no points are cut off, use breaks as is
  # if clip==FALSE and no points are cut off, use breaks as is

  # X breaks
  xclip_min <- any(naf(data$xorg < xylimits$xlim["min"]))
  xclip_max <- any(naf(data$xorg > xylimits$xlim["max"]))

  xbreaks <- scales::pretty_breaks(n = nbreaks_x)(xylimits$xlim, n = nbreaks_x)
  if (clip & xclip_min){
    xylimits$xlim["min"] <- min(xbreaks)
    xclip_min <- any(data$xorg < xylimits$xlim["min"])
  }
  if (clip & xclip_max){
    xylimits$xlim["max"] <- max(xbreaks)
    xclip_max <- any(data$xorg > xylimits$xlim["max"])
  }


  xylimits$xlim <- xylimits$xlim + c(-diff(xylimits$xlim), diff(xylimits$xlim)) * c(!xclip_min, !xclip_max)*0.02
  # if (xclip_min){ xbreaks[xbreaks < xylimits$xlim[1]] <- xylimits$xlim[1] } # adjust breaks?
  # if (xclip_max){ xbreaks[xbreaks > xylimits$xlim[2]] <- xylimits$xlim[2] }

  # names(xbreaks) <- paste0(ifelse(xbreaks < 0, "-", ""), abs(xbreaks))
  names(xbreaks) <- as.character(xbreaks)
  if (xclip_min){ names(xbreaks)[1] <- paste0("<", xbreaks[1]) }
  if (xclip_max){ names(xbreaks)[length(xbreaks)] <- paste0(">", xbreaks[length(xbreaks)]) }



  # Y breaks
  yclip_min <- any(naf(data$yorg < xylimits$ylim["min"]))
  yclip_max <- any(naf(data$yorg > xylimits$ylim["max"]))

  ybreaks <- scales::pretty_breaks(n = nbreaks_y)(xylimits$ylim, n = nbreaks_y)

  if (clip & yclip_min){
    xylimits$ylim["min"] <- min(ybreaks)
    yclip_min <- any(data$xorg < xylimits$ylim["min"])
  }
  if (clip & yclip_max){
    xylimits$ylim["max"] <- max(ybreaks)
    yclip_max <- any(data$xorg > xylimits$ylim["max"])
  }


  xylimits$ylim <- xylimits$ylim + c(-diff(xylimits$ylim), diff(xylimits$ylim)) * c(!yclip_min & !at_zero, !yclip_max)*c(0.01, 0.05)
  # if (yclip_min){ ybreaks[ybreaks < xylimits$ylim[1]] <- xylimits$ylim[1] }
  # if (yclip_max){ ybreaks[ybreaks > xylimits$ylim[2]] <- xylimits$ylim[2] }

  # names(ybreaks) <- paste0(ifelse(ybreaks < 0, "-", ""), abs(ybreaks))
  names(ybreaks) <- as.character(ybreaks)
  if (yclip_min){ names(ybreaks)[1] <- paste0("<", ybreaks[1]) }
  if (yclip_max){ names(ybreaks)[length(ybreaks)] <- paste0(">", ybreaks[length(ybreaks)]) }


  data$x[ data$x < xylimits$xlim["min"] ] <- xylimits$xlim["min"]
  data$x[ data$x > xylimits$xlim["max"] ] <- xylimits$xlim["max"]
  data$y[ data$y < xylimits$ylim["min"] ] <- xylimits$ylim["min"]
  data$y[ data$y > xylimits$ylim["max"] ] <- xylimits$ylim["max"]


  ### GGPLOT ###

  data <- data[order(data$score, decreasing = FALSE),]

  gg <- data %>% ggplot2::ggplot(ggplot2::aes(x = x, y = y, label = label, color = !!color, shape = !!shape, ...))

  if (show_grid){
    panel_grid <- ggplot2::element_line(size = lwd, color = rgb(0.9,0.9,0.9))
  } else {
    panel_grid <- ggplot2::element_blank()
  }

  gg %<>% + ggplot2::theme_bw(base_size = 20)
  gg %<>% + ggplot2::theme(text = ggplot2::element_text(color = "black", size = lab_size),
                  rect = ggplot2::element_rect(color = "black", size = lwd),
                  line = ggplot2::element_line(size = lwd),
                  legend.text = ggplot2::element_text(color = "black", size = leg_size),
                  legend.title = ggplot2::element_text(color = "black", size = leg_size),
                  panel.grid.minor = ggplot2::element_blank(),
                  panel.grid.major = panel_grid,
                  panel.border = ggplot2::element_rect(colour = "black", fill = NA, size = ifelse(clip_frame == "on", lwd*2, lwd)),
                  strip.background = ggplot2::element_blank(),
                  strip.text = ggplot2::element_text(color = "black", size = title_size),
                  axis.ticks = ggplot2::element_line(color = "black", size = lwd),
                  axis.line = ggplot2::element_blank(),
                  plot.margin = ggplot2::unit(c(1,1,1,1), "cm"),
                  plot.title = ggplot2::element_text(size = title_size, hjust = 0.5, lineheight = 1.5),
                  axis.title = ggplot2::element_text(size = axis_size, face = "plain"),
                  axis.text = ggplot2::element_text(size = axis_size, color = "black"))

  if (!is.null(ptres)){ gg %<>% + ggplot2::geom_hline(yintercept = -log10(ptres), linetype = "dashed", color = rgb(0.3,0.3,0.3)) }

  # points
  if (scale_size == FALSE){
    gg %<>% + ggplot2::geom_point(size = point_size, alpha = 0.8, stroke = stroke)
  } else {
    gg %<>% + ggplot2::geom_point(aes(size = score), alpha = 0.8, stroke = stroke)
    gg %<>% + ggplot2::scale_size_continuous(range = c(point_size/5, point_size*2), guide = "none")
  }

  if (color_user_def){
    color_vals <-  setNames(c(color_up, color_down, color_nonsig), c(label_up, label_down, label_nonsig))
    if (show_nonsig){
      color_breaks <- c(label_up, label_nonsig, label_down)
    } else {
      color_breaks <- c(label_up, label_down)
    }
    gg %<>% + ggplot2::scale_colour_manual(values = color_vals, breaks = color_breaks)
    gg %<>% + ggplot2::labs(title = title, y = paste0("-log10 ", rlang::as_name(y)), x = rlang::as_name(x), color = NULL)
  } else {
    gg %<>% + ggplot2::labs(title = title, y = paste0("-log10 ", rlang::as_name(y)), x = rlang::as_name(x))
  }

  gg %<>% + ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0,0)),
                               limits = xylimits$xlim,
                               breaks = xbreaks,
                               labels = names(xbreaks))

  gg %<>% + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,0)),
                               limits = xylimits$ylim,
                               breaks = ybreaks,
                               labels = names(ybreaks))

  # point labels
  if (is.null(attract)) attract <- sqrt(repel)
  gg %<>% + ggrepel::geom_text_repel(data = subset(data, do_label == TRUE),
                                     fontface = labface,
                            size = lab_size/ggplot2:::.pt,
                            seed = seed,
                            xlim = xylimits$xlim - c(-diff(xylimits$xlim), diff(xylimits$xlim))*0.18,
                            ylim = xylimits$ylim - c(-diff(xylimits$ylim)*0.3, diff(xylimits$ylim)*0.02),
                            force = repel, force_pull = attract,  max.overlaps = max_overlaps,
                            point.padding = 0.35, box.padding = box.padding, max.time = 30, max.iter = 10^6,
                            min.segment.length = 0, vjust = 0, color = rgb(0.0,0.0,0.0), segment.alpha = 0.6)

  gg %<>% + ggplot2::guides(size = "none", label = "none", color = ggplot2::guide_legend(override.aes = list(size = leg_key_size)))
  gg %<>% + ggplot2::coord_cartesian(clip = clip_frame)

  return(gg)
}





getLimits <- function(x, clip = TRUE, expand = 0.1, negative = TRUE){

  x <- x[!is.na(x)]
  x <- x + x*expand

  if (clip == TRUE){
    h <- hist(x, plot = FALSE, breaks = 30)
    xd <- (h$counts > 3) | (rev(cumsum(rev(h$counts))) > 8)
    xmin <- h$breaks[which(xd)[1]]
    xmax <- rev(h$breaks)[which(rev(xd))[1]]

  } else {
    xmin <- NA
    xmax <- NA
  }

  # if (is.na(xmin)){ xmin <- floor(min(x)*10^-floor(-log10(abs(min(x)))))/10^-floor(-log10(abs(min(x)))) }
  # if (is.na(xmin)){xmin <- 0}
  # if (xmin > 0){xmin <- 0}
  # if (is.na(xmax)){ xmax <- ceiling(max(x)*10^-floor(-log10(abs(max(x)))))/10^-floor(-log10(abs(max(x)))) }

  if (is.na(xmax)) xmax <- max(x) %>% roundup(., roundup(-log10(abs(.)))) # upper
  if (is.na(xmin)) xmin <- min(x) %>% rounddown(., roundup(-log10(abs(.)))) # lower

  if (is.na(xmin)){xmin <- -0.1 * xmax}
  if (is.na(xmax)){xmax <- -0.1 * xmin}

  if (is.na(xmax) & is.na(xmin)){
    xmin <- -1
    xmax <- 1
  }


  res <- c("min" = xmin, "max" = xmax)
  if (negative == FALSE) res[res < 0] <- 0

  res
}





















