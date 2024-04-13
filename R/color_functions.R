



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

  if (class(colors) == "list"){
    lapply(colors, prismatic::color)
  } else {
    prismatic::color(colors)
  }
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










#' Get color gradients
#'
#' @param colors vector of colors
#' @param breaks vector of breaks (required for 'colorRamp2')
#' @param type generate gradient as ggplot2, colorRamp2 or none
#' @param name use a pre-defined gradient
#'
#' @return
#' @export
#'
#' @examples
#' gradient(colors = c("-3 "= "blue", "0" = "white", "3" = "orangered"))
#' gradient(name = "bluered3", type = "colorramp2")
#' ggplot(iris, aes(Sepal.Length, Petal.Length, color = Sepal.Width)) + geom_point() + gradient(c("red", "green", "blue"), breaks = c(0,1,5))
#' ggplot(iris, aes(Sepal.Length, Petal.Length, color = Sepal.Width)) + geom_point() + gradient(c("red", "green", "blue"))
#' ComplexHeatmap::Heatmap(iris[,-5], col = gradient(name = "blueyellow3", type = "colorramp2", breaks = c(0,4,8)))
gradient <- function(colors, breaks = NULL, type = c("ggplot2", "colorRamp2"), name = NULL){

  if (!is.null(name)){
    if (name == "bluered3") colors <- c("blue", "white", "red")
    if (name == "blueyellow3") colors <- c("#146bc7", "white", "#eb9d0e")
  }

  if (is.null(breaks) & !is.null(names(colors))){
    breaks <- as.numeric(names(colors))
  }

  # convert to output type
  type <- type[1]

  if (tolower(type) == "ggplot2"){
    if (is.null(breaks)){
      ggplot2::scale_color_gradientn(colours = colors)
    } else {
      ggplot2::scale_color_gradientn(colours = colors, values = rangescale(breaks), limits = c(min(breaks), max(breaks)))
    }

  } else if (tolower(type) == "colorramp2"){
    if (is.null(breaks)){
      warning("Warning: Automatically setting breaks for 'colorRamp2' color scale")
      breaks <- seq(colors) - mean(seq(colors)) * 2
    }
    circlize::colorRamp2(breaks = as.numeric(breaks), colors = colors)

  } else {
    colors

  }
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






#' Select colors using the 'paletteer' package
#'
#' @param length Integer/vector of desired length(s)
#' @param continuous TRUE
#' @param discrete TRUE
#'
#' @return
#' @export
#'
#' @examples
#' getPaletteer(length = 2)$colors %>% showColors()
getPaletteer <- function(length = NULL, continuous = TRUE, discrete = TRUE){

  stopifnot(requireNamespace("paletteer"))

  cont <- paletteer::palettes_c_names
  disc <- paletteer::palettes_d_names

  cont$class <- "continuous"
  cont$length <- NA
  disc$class <- "discrete"
  disc$novelty <- NULL
  df <- rbind(disc, cont[,colnames(disc)])

  if (!is.null(length)){
    df <- df[df$length %in% length,, drop = FALSE]
  }

  if (continuous == FALSE){
    df <- df[df$class != "continuous",, drop = FALSE]
  }

  if (discrete == FALSE){
    df <- df[df$class != "discrete",, drop = FALSE]
  }

  colors <- unique(df$package) %L>% function(pkg){
    tmp <- df[df$package == pkg,]
    tmp$palette %L>% function(pal){
      paletteer::palettes_d[[pkg]][[pal]]
    }
  }
  df$colors <- unlist(colors, recursive = FALSE)

  colorcat("Use 'showColors(df$colors)' to display colors.", col = rgb(0.3,0.3,1))
  df
}














#' Basic ggplot2 theme
#'
#' @param fontsize Base font size
#' @param fontfamily Base font family
#' @param lwd Base linewidth
#' @param color Base color
#' @param grid.color Grid color (no grid if NULL)
#' @param title.face Title face (plain/bold/...)
#' @param title.cex Title size factor
#' @param axis.title.cex Axis title size factor
#' @param axis.text.cex Axis text size factor
#' @param grid.cex Grid linewidth size factor
#' @param ... other parameters passed to theme(...)
#' @inheritParams ggplot2::theme_bw
#'
#' @export
#'
#' @examples
#' ggplot(iris, aes(Sepal.Length, Petal.Length)) + geom_point() + ggtitle("Iris") + theme_basic()
#' ggplot(iris, aes(Sepal.Length, Petal.Length, color = Species)) + geom_point() + ggtitle("Iris") + theme_basic(grid.color = "grey70")
theme_basic <- function(fontsize = 18, fontfamily = "", lwd = NULL, color = "black", grid.color = NULL,
                        title.face = "bold", title.cex = 1, axis.title.cex = 1, axis.text.cex = 0.8, grid.cex = 0.7, ...){

  if (is.null(lwd)) lwd <- fontsize/22

  th1 <- ggplot2::theme_bw(base_size = fontsize,
                           base_family = fontfamily,
                           base_line_size = lwd,
                           base_rect_size = lwd)

  th2 <- ggplot2::theme(line = ggplot2::element_line(colour = color, linewidth = lwd),
                        rect = ggplot2::element_rect(colour = color, fill = NA, linewidth = lwd),
                        text = ggplot2::element_text(colour = color, size = fontsize, family = fontfamily),
                        title = ggplot2::element_text(colour = color, size = fontsize, family = fontfamily),
                        axis.text = ggplot2::element_text(colour = color, size = ggplot2::rel(axis.text.cex)),
                        axis.title = ggplot2::element_text(colour = color, size = ggplot2::rel(axis.title.cex)),
                        axis.ticks = ggplot2::element_line(colour = color),
                        axis.line = ggplot2::element_line(colour = color, lineend = "square"),
                        legend.text = ggplot2::element_text(colour = color),
                        legend.title = ggplot2::element_text(colour = color, hjust = 0),
                        panel.border = ggplot2::element_blank(),
                        panel.grid = ggplot2::element_blank(),
                        plot.title = ggplot2::element_text(colour = color, face = title.face, hjust = 0.5, size = ggplot2::rel(title.cex)),
                        strip.background = ggplot2::element_blank(),
                        strip.text = ggplot2::element_text(colour = color, size = fontsize, family = fontfamily),
                        strip.clip = "off",
                        ...)

  if (!is.null(grid.color)){
    th2 <- ggplot2::`%+replace%`(th2, ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = ggplot2::rel(grid.cex), colour = grid.color)))
  }

  th <- ggplot2::`%+replace%`(th1, th2)
  th

}







#' Ggplot2 theme for checking plot layout and margins
#'
#' @param theme_orig
#'
#' @return
#' @export
#'
#' @examples
theme_dev <- function(theme_orig){
  th2 <- theme(panel.background = element_rect(fill = rgb(0.5, 0.5, 0.5)),
               plot.background = element_rect(fill = "darkblue"))
  ggplot2::`%+replace%`(theme_orig, th2)
}






