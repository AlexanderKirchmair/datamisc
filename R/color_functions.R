



#' Interactive color picker
#'
#' @export
#'
pickColors <- function(...){
  colourpicker::colourWidget(...)
}




#' Display colors in console
#'
#' @param colors input colors
#'
#' @return object of class colors
#' @export
#'
#' @examples
showColors <- function(colors){
  if (missing(colors)) colors <- c("red", "green", "blue", sample(grDevices::colors(), 6))
  prismatic::color(colors)
}




#' Print in color
#'
#' @param str
#' @param col
#' @param add
#'
#' @return
#' @export
#'
#' @examples
colorcat <- function(str = "text", col = rgb(1,1,1), add = "\n"){
  cat( crayon::make_style(col)(paste0(str, add)) )
}





#' Generate color mappings for discrete data columns in a dataframe
#'
#' @param dataframe
#'
#' @return
#' @export
#'
#' @examples
getColors <- function(dataframe){

  if (class(dataframe) != "data.frame") dataframe <- data.frame(dataframe)

  if (is.null(dataframe)) return(NULL)
  if (ncol(dataframe) == 0) return(NULL)

  colors <- list()


  ### Discrete
  dataframe.discrete <- dataframe[, sapply(dataframe, is.factor) | sapply(dataframe, is.character), drop=FALSE]
  n <- sapply(dataframe.discrete, function(x) length(unique(x)))

  colors <- genPalettes(n = NULL, length_each = n)

  colors <- lapply(seq_along(colors), function(i){
    tmp <- colors[[i]]
    names(tmp) <- unique(dataframe.discrete[[i]])
    tmp
  })

  names(colors) <- colnames(dataframe.discrete)

  colors
}



#' Generate distinct color palettes
#'
#' @param n Number of palettes
#' @param length_each Length of each palette (either a single digit or a vector of length n)
#' @param dist Distance between palettes
#' @param saturation
#' @param lightness
#' @param cvd
#' @param cvd_severity
#' @param seed
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' genPalettes(n = 10) %L>% showColors()
genPalettes <- function(n = 1, length_each = 3, dist = 0.75, saturation = c(0.5, 0.7), lightness = c(0.5, 0.7), cvd = "protan", cvd_severity = 0.2, seed = 123, ...){

  if (is.null(n)){
    n <- length(length_each)
  }

  if (length(length_each) == 1 & n != 1){
    length_each <- rep(length_each, n)
  }

  if (length(length_each) != n){
    stop("Please provide the number of colors to generate in each palette!")
  }

  # generate colorspaces for each palette
  rel <- (0:n)/(n)
  drel <- mean(diff(rel))
  sep <- 1 - dist
  init <- scales::rescale(rel, from = c(min(rel)-drel*sep, max(rel)+drel*sep), to = c(1, 359))
  d <- mean(diff(init))

  set.seed(seed)
  init <- init[-length(init)] + d * runif(1)
  spaces <- lapply(init, function(tmp){
    list(h = c(tmp - d * sep, tmp + d * sep), s = saturation, l = lightness)
  })

  # generate palette colors
  pals <- lapply(1:n, function(i){
    qualpalr::qualpal(n = length_each[i], colorspace = spaces[[i]], cvd_severity = cvd_severity, cvd = cvd, ...)$hex
  })

  # sort by approximate brightness
  pals <- lapply(pals, function(pal){
    ix <- order(grDevices::rgb2hsv(grDevices::col2rgb(pal))["v",], decreasing = TRUE)
    pal[ix]
  })

  names(pals) <- names(length_each)
  pals
}









#' Basic ggplot2 theme
#'
#' @param base_size
#' @param base_family
#' @param base_line_size
#' @param base_rect_size
#' @param base_color
#' @param grid
#' @param ... other parameters passed to theme(...)
#' @inheritParams ggplot2::theme_bw
#'
#' @export
#'
#' @examples
theme_basic <- function(base_size = 18, base_family = "", base_line_size = base_size/22, base_rect_size = base_size/22, base_color = "black", grid = FALSE, ...){

  th1 <- ggplot2::theme_bw(base_size = base_size,
                          base_family = base_family,
                          base_line_size = base_line_size,
                          base_rect_size = base_rect_size)

  th2 <- ggplot2::theme(line =  ggplot2::element_line(colour = base_color, size = base_line_size),
                   rect =  ggplot2::element_rect(colour = base_color, fill = NA, size = base_rect_size),
                   text =  ggplot2::element_text(colour = base_color, size = base_size, family = base_family),
                   title =  ggplot2::element_text(colour = base_color),
                   axis.text =  ggplot2::element_text(colour = base_color, size =  ggplot2::rel(0.75)),
                   axis.ticks =  ggplot2::element_line(colour = base_color),
                   axis.line =  ggplot2::element_line(colour = base_color),
                   legend.text =  ggplot2::element_text(colour = base_color),
                   legend.title =  ggplot2::element_text(colour = base_color),
                   legend.title.align = 0,
                   panel.border =  ggplot2::element_blank(),
                   panel.grid =  ggplot2::element_blank(),
                   plot.title =  ggplot2::element_text(colour = base_color),
                   strip.background =  ggplot2::element_blank(),
                   strip.text =  ggplot2::element_text(colour = base_color),
                   ...)

  th <- ggplot2::`%+replace%`(th1, th2)

  th

}














