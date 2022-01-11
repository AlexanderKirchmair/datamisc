


# NOTE
# colors should be stored in a list as following:
# colors$factor = c(lev1 = "blue", lev2 = "red)
# colors$continuous = c("min" = "blue", "mid" = "white", "max" = "red")



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



#' Generate color mappings for discrete data columns in a dataframe
#'
#' @param dataframe
#'
#' @return
#' @export
#'
#' @examples
getColors <- function(dataframe){

  ### Function to generate colors for columns in dataframe

  if (class(dataframe) != "data.frame") dataframe <- data.frame(dataframe)

  stopifnot(require("paletteer"))

  if (is.null(dataframe)) return(NULL)
  if (ncol(dataframe) == 0) return(NULL)

  colors <- list()


  ### Discrete
  dataframe.discrete <- dataframe[,sapply(dataframe, is.factor) | sapply(dataframe, is.character),drop=FALSE]

  pals <- subset(palettes_d_names, type == "qualitative")
  pals <- pals[order(pals$length),]
  pals$used <- FALSE

  for (i in seq_along(colnames(dataframe.discrete))){

    elem <- unique(dataframe.discrete[[i]])
    usepal <- which(pals$length >= length(elem) & !pals$used)[1]

    if (!is.na(usepal)){
      pals$used[usepal] <- TRUE
      tmpcol <- palettes_d[[pals$package[usepal]]][[pals$palette[usepal]]]
      tmpcol <- tmpcol[1:length(elem)]
      colors[[colnames(dataframe.discrete)[i]]] <- setNames(tmpcol, elem)
    }
  }

  return(colors)
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

  th <- ggplot2::theme_bw(base_size = base_size,
                          base_family = base_family,
                          base_line_size = base_line_size,
                          base_rect_size = base_rect_size)

  th <- th %+replace%
    ggplot2::theme(line = element_line(colour = base_color, size = base_line_size),
                   rect = element_rect(colour = base_color, fill = NA, size = base_rect_size),
                   text = element_text(colour = base_color, size = base_size, family = base_family),
                   title = element_text(colour = base_color),
                   axis.text = element_text(colour = base_color, size = rel(0.75)),
                   axis.ticks = element_line(colour = base_color),
                   axis.line = element_line(colour = base_color),
                   legend.text = element_text(colour = base_color),
                   legend.title = element_text(colour = base_color),
                   legend.title.align = 0,
                   panel.border = element_blank(),
                   panel.grid = element_blank(),
                   plot.title = element_text(colour = base_color),
                   strip.background = element_blank(),
                   strip.text = element_text(colour = base_color),
                   ...)


  th

}

#
# pal <- qualpalr::autopal(10, colorspace = "pretty_dark")
# showColors(pal$hex)














