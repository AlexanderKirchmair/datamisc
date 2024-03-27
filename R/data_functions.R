

### Basic utils ------


#' Return the first n rows and columns of an object
#'
#' @param data
#' @param nrows
#' @param ncols
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
head2 <- function(data, nrows = 15, ncols = 10, ...){
  rows <- ifelse(nrow(data) < nrows, nrow(data), nrows)
  columns <- ifelse(ncol(data) < ncols, ncol(data), ncols)
  data[1:rows, 1:columns, drop = FALSE]
}



#' Return the middle n rows and columns of an object
#'
#' @param data
#' @param nrows
#' @param ncols
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' body(mtcars)
body <- function(data, nrows = 15, ncols = 10, ...){

  nr <- nrow(data)
  nc <- ncol(data)

  rr <- round(mean(1:nr)) + c(-ceiling(nrows/2), floor(nrows/2)-1)
  if (min(rr) < 1) rr <- rr - min(rr) + 1
  if (max(rr) > nr) rr[2] <- nr

  rc <- round(mean(1:nc)) + c(-ceiling(ncols/2), floor(ncols/2)-1)
  if (min(rc) < 1) rc <- rc - min(rc) + 1
  if (max(rc) > nc) rc[2] <- nc

  data[rr[1]:rr[2], rc[1]:rc[2], drop = FALSE]

}


#' Limit number of characters in a string
#'
#' @param x character string
#' @param maxchar max. number of characters
#' @param add string added to truncated objects
#' @param add_incl include added string in maxchar
#'
#' @return
#' @export
#'
#' @examples
cutstr <- function(x, maxchar = 25, add = "...", add_incl = TRUE){
  ix <- nchar(x) > maxchar
  x[ix] <- substr(x[ix], 1, ifelse(add_incl, maxchar - nchar(add), maxchar))
  x[ix] <- paste0(x[ix], add)
  x
}


#' File extension of path
#'
#' @param path
#'
#' @export
#'
#' @examples
#' baseext(list.files())
#'
baseext <- function(path, ...){
  path <- basename(path)
  ext <- gsub(x = path, pattern = ".*\\.", replacement = "")
  ext[!grepl(pattern = ".", x = path, fixed = TRUE)] <- ""
  ext
}




#' File path
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fp <- function(...){
  args <- list(...)
  args <- lapply(args, paste0, collapse = "")
  do.call(file.path, args)
}



#' Memory size of workspace objects
#'
#' @param x
#' @param units
#'
#' @return
#' @export
#'
#' @examples size(ls())
size <- function(x, units = NULL){

  if (is.character(x)) x <- lapply(x, get)

  s <- object.size(x)
  o <- log10(s)

  if (is.null(units)){
    if (o > 2) units <- "Kb"
    if (o > 5) units <- "Mb"
    if (o > 8) units <- "Gb"
  }

  format(s, units = units)
}


#' Get class of a quoted object
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
quo_class <- function(x, ...){
  x_class <- NA
  try({x_class <- class(rlang::eval_tidy(x, ...))}, silent = TRUE)
  x_class
}


#' Divide vector into n groups
#'
#' @param x
#' @param n
#'
#' @return
#' @export
#'
#' @examples
#' partition(LETTERS[1:10], 3)
partition <- function(x, n = 3){
  stopifnot(length(x) >= n)
  if (n == 1) return(rep(1, length(x)))
  y1 <- rep(1:(n-1), each = floor(length(x)/n))
  y2 <- rep(n, length(x) - length(y1))
  c(y1, y2)
}


#' Pattern Matching and Replacement
#'
#' @param pattern Multiple patterns
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
lgrep <- function(pattern, x, ...){
  ix <- sapply(pattern, grep, x = x, ...)
  sort(unique(unlist(ix)))
}


#' Colored printing to compare two dataframes
#'
#' @param df1
#' @param df2
#' @param ncol
#' @param signif
#' @param sep
#' @param marker
#' @param color
#'
#' @return
#' @export
#'
#' @examples
compareDF <- function(df1, df2, ncol = 10, signif = 3, sep = "/", marker = "*", color = "red"){

  m1 <- data.matrix(df1)
  m2 <- data.matrix(df2)

  stopifnot(all(dim(m1) == dim(m2)))

  diff <- ifelse(naf(m1 == m2), "", marker)
  diff[is.na(m1) | is.na(m2)] <- marker
  diff[is.na(m1) & is.na(m2)] <- ""

  com <- paste0(nar(as.character(unlist(signif(m1, signif))), "NA"), sep, nar(as.character(unlist(signif(m2, signif))), "NA"))
  com <- paste0(com, diff)

  ix <- !grepl(marker, x = com, fixed = TRUE)
  com[ix] <- m1[ix]

  out <-  matrix(com, nrow = nrow(m1), dimnames = dimnames(m1))

  if (!is.null(ncol)) out <- out[,1:ncol]

  cdf <- data.frame(out)
  cdf <- colorDF::df_search(cdf, pattern = "\\*")

  attributes(cdf)$.style$id <- "universal"
  attributes(cdf)$.style$sep <- NULL
  attributes(cdf)$.style$col.names$bg <- NULL
  attributes(cdf)$.style$col.names$fg <- "grey60"
  attributes(cdf)$.style$row.names$fg <- "grey60"
  attributes(cdf)$.style$col.names$align <- "center"
  attributes(cdf)$.style$interleave$bg <- NULL
  attributes(cdf)$.style$row.names$decoration <- "NULL"
  attributes(cdf)$.style$type.styles$match$align <- "right"
  attributes(cdf)$.style$type.styles$match$fg_match <- color

  colorDF::print_colorDF(cdf)
}





#' Column to rownames
#'
#' @param data
#' @param col
#' @param sep
#'
#' @return
#' @export
#'
#' @examples
col2rownames <- function(data, col = id, sep = "_"){

  col <- rlang::enquo(col)

  if (!rlang::as_name(col) %in% colnames(data)){
    warning(paste0("Warning: Column ", rlang::as_name(col), " not found!"))
    return(data)
  }

  names <- dplyr::select(data, !!col)
  names <- apply(names, 1, paste0, collapse = sep)
  rownames(data) <- names
  data <- dplyr::select(data, -!!col)
  data
}




#' Rownames to column
#'
#' @param data
#' @param col
#' @param keep
#'
#' @return
#' @export
#'
#' @examples
rownames2col <- function(data, col = id, keep = FALSE){

  col <- rlang::enquo(col)

  names <- rownames(data)
  if (is.null(names)){
    warning("Warning: No rownames found.")
    return(data)
  }
  if (keep == FALSE) rownames(data) <- NULL

  i <- rlang::as_name(col)
  if (i %in% colnames(data)){
    data[,i] <- names
  } else {
    if (class(names) %in% unique(sapply(data, class))){
      data <- cbind(names, data)
    } else {
      data <- data.frame(names, data)
    }

    colnames(data)[1] <- i
  }

  data
}







### NA handling ------

#' Set NA values to FALSE
#' @return
#' @export
naf <- function(data, ...){
  data[is.na(data)] <- FALSE
  data
}


#' Set NA values to TRUE
#' @return
#' @export
nat <- function(data, ...){
  data[is.na(data)] <- TRUE
  data
}

#' Set NA values to any other value
#' @return
#' @export
nar <- function (data, replace = 0, ...){
  data[is.na(data)] <- replace
  data
}





#' Skip rows/columns containing NA values
#'
#' Function for fast subsetting of data (see also 'na.omit()').
#'
#'
#' @param data Dataset
#' @param skip Skip rows (columns) containing NA values.
#'
#' @return
#' @export
na_skip <- function(data, skip = "rows"){
  if (tolower(skip) %in% c("row", "rows")) i <- 1
  else if (tolower(skip) %in% c("col", "cols", "columns")) i <- 2
  else stop("Error: skip must be 'rows' or 'columns'.")

  ix <- !apply(is.na(data), i, any)

  if (i == 1) data <- data[ix,,drop = FALSE] else data <- data[,ix,drop = FALSE]

  data
}




### Random data generation ------


#' Generate a matrix of random numbers
#'
#' @description
#' Returns a matrix of size nrow x ncol from random numbers generated by FUN.
#'
#' @param nrow number of rows
#' @param ncol number of columns
#' @param FUN random number generator
#' @param ... parameters other than n passed to FUN
#'
#' @return
#' @export
#'
#' @examples
#' @seealso \code{\link[stats]{runif}}, \code{\link[stats]{rnorm}}, \code{\link[stats]{rexp}}
rmat <- function(nrow = 3, ncol = 5, FUN = runif, ...){
  n <- nrow * ncol
  v <- FUN(n = n, ...)
  matrix(data = v, nrow = nrow, ncol = ncol,
         dimnames = list(paste0("r", 1:nrow), paste0("c", 1:ncol)))
}


#' Generate a random list
#'
#' @param length length of the list
#' @param items (variable) number of items per list entry
#' @param space set of items to sample
#'
#' @return list
#' @export
#'
#' @examples
#' rlist()
#' rlist(length = 3, items = 2:4, space = 1:10)
rlist <- function(length = 5, items = 1:3, space = LETTERS){
  if (length(space) < max(items)) stop("Error: Requested number of items cannot exceed number of items in space!")
  res <- lapply(1:length, function(x) sample(space, size = sample(c(items, items), size = 1)) )
  setNames(res, paste0(seq(res), "_", sapply(res, function(x) paste0(x, collapse = "") ) ))
}





#' Generate a random data.frame
#'
#' @param nrow
#' @param ncol
#'
#' @return
#' @export
#'
#' @examples
#' rdataframe()
rdataframe <- function(nrow = 5, ncol = 5){

  df <- as.data.frame(lapply(1:ncol, function(i){
    if (runif(1) > 0.5){
      sample(sample(LETTERS, max(2, round(nrow/4))), size = nrow, replace = TRUE)
    } else {
      sample(as.character(1:(round(nrow/2))), size = nrow, replace = TRUE)
    }
  }))

  colnames(df) <- sample(LETTERS, ncol(df))
  df
}







### Duplicate data manipulation ------


#' Get unique sets (irrespective of item order)
#'
#' @param sets list of sets
#' @param sep collapse set names
#'
#' @return
#' @export
#'
#' @examples
uniqueSets <- function(sets, sep = "&"){

  nonames <- is.null(names(sets))
  if (nonames) names(sets) <- seq(sets)

  which_dups <- sapply(sets, function(tmp1) sapply(sets, function(tmp2) {setequal(tmp1, tmp2) } ))
  diag(which_dups) <- FALSE

  dup_names <- sapply(colnames(which_dups), function(tmp1) rownames(which_dups)[which_dups[,tmp1]] )
  dedup_names <- sapply(names(dup_names), function(tmp){
    if (length(dup_names[[tmp]]) > 1){
      paste(sort(c(tmp, dup_names[[tmp]])), collapse = sep)
    } else {
      tmp
    }
  })

  names(sets) <- dedup_names
  sets <- sets[unique(dedup_names)]

  if (nonames) sets <- setNames(sets, NULL)
  sets <- lapply(sets, sort)
  return(sets)
}


#' Rename duplicated strings
#'
#' @param x character vector
#' @param sep separator
#' @param index indices added to duplicated elements
#' @param ... arguments passed to 'duplicated'
#'
#' @return
#' @export
#'
#' @examples
#' replicate(3, paste(LETTERS[1:3], collapse = "")) %>% dedupl()
#' replicate(3, paste(LETTERS[1:3], collapse = "")) %>% dedupl(sep = "_", index = letters)
#'
dedupl <- function(x, sep = ".", index = NULL, ...){

  xorg <- x
  ix <- duplicated(x, ... = )
  i <- 1

  if (is.null(index)){
    while (any(ix)){
      x[ix] <- paste0(xorg[ix], sep, i)
      ix <- duplicated(x, ...)
      i <- i + 1
    }
  } else {
    while (any(ix)){
      x[ix] <- paste0(xorg[ix], sep, index[i])
      ix <- duplicated(x, ...)
      i <- i + 1
    }
  }

  stopifnot(!any(duplicated(x, ...)))
  x
}




### Numeric data manipulations ------


#' Scale rows (columns) of a matrix
#'
#' @description
#' Apply FUN ('scale' by default) to the rows or columns of a numeric matrix.
#'
#' @param data Matrix or data.frame
#' @param rows Scale rows (TRUE/FALSE)
#' @param cols Scale columns (TRUE/FALSE)
#' @param FUN Function used for scaling
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' matScale(rmat(5, 5), rows = TRUE)
matScale <- function(data, rows = FALSE, cols = FALSE, FUN = scale, ...){

  data.org <- data
  ix <- sapply(as.data.frame(data), is.numeric)
  data <- data[,ix, drop = FALSE]
  names.org <- dimnames(data)

  if (rows == TRUE & cols == TRUE) stop("Error: Do not scale rows and columns at once!")

  if (rows == TRUE){
    data <- t(apply(data, 1, function(tmp) as.numeric(FUN(tmp, ...))))
  }

  if (cols == TRUE){
    data <- apply(data, 2, function(tmp) as.numeric(FUN(tmp, ...)))
  }

  dimnames(data) <- names.org
  data <- cbind(data, data.org[,!ix, drop = FALSE])
  if (!is.null(colnames(data.org))) data <- data[,colnames(data.org)]

  stopifnot( all.equal(dim(data.org), dim(data)) )
  return(data)
}


#' Replace Inf values for plotting
#'
#' @param data
#' @param increase
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
subInf <- function(data, increase = 0.2, ...){

  data.noinf <- data
  data.noinf[is.infinite(data)] <- NA
  maxval <- roundup(max(abs(data.noinf), na.rm = TRUE)*increase)

  data[data == Inf] <- maxval
  if (any(data < 0, na.rm = TRUE)) data[data == -Inf] <- -maxval

  data
}


#' Rounding of numbers
#'
#' @description
#' Round up
#'
#' @param x numeric vector
#' @param digits number of decimal places
#' @param ... other arguments passed to ceiling/floor
#'
#' @export
#'
roundup <- function(x, digits = 0, ...){
  ceiling(x * 10^digits, ...) / 10^digits
}


#' Rounding of numbers
#'
#' @description
#' Round down
#'
#' @inheritParam
#'
#' @export
#'
rounddown <- function(x, digits = 0, ...){
  floor(x * 10^digits, ...) / 10^digits
}




#' Adjust P-values for Multiple Comparisons
#'
#' @description
#' Adjust a vector or matrix of p-values derived from the same data for multiple testing.
#'
#' @param p vector or matrix of p-values
#' @param ...
#'
#' @inheritParam stats::p.adjust
#'
#' @return
#' @export
#'
padjust <- function(p, method = "fdr", ...){
  porig <- p
  is.mat <- !is.null(dim(porig))
  if (is.mat) p <- as.vector(data.matrix(p))
  padj <- stats::p.adjust(p, method = method, ...)
  if (is.mat) padj <- matrix(padj, nrow = nrow(porig), dimnames = dimnames(porig))
  if ("data.frame" %in% class(porig)) padj <- as.data.frame(padj)
  padj
}



### Dataframe processing ------


#' Cbind multiple matrices by shared rownames
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
cjoin <- function(...){

  quos_args <- rlang::quos(...)
  dfnames <- sapply(quos_args, rlang::quo_name)

  dfnames1 <- names(list(...))
  if(!is.null(dfnames1)) dfnames[nchar(dfnames1)>0] <- dfnames1[nchar(dfnames1)>0]

  dflist <- lapply(quos_args, rlang::eval_tidy)
  dfs <- lapply(dflist, as.data.frame)

  for (i in which(sapply(dfs, ncol) == 1)){ colnames(dfs[[i]]) <- dfnames[i] }

  allids <- lapply(dfs, rownames)
  ids <- Reduce(f = intersect, x = allids)

  cols <- unlist(lapply(dfs, colnames))
  if (any(duplicated(cols))) cols <- paste0(cols, ".", rep(dfnames, sapply(dfs, ncol)))

  subdfs <- lapply(dfs, function(tmp) tmp[ids,,drop = FALSE] )
  res <- Reduce(f = cbind, x = subdfs)
  colnames(res) <- cols
  res
}


#' Rbind multiple matrices by shared colnames
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
rjoin <- function(...){

  quos_args <- rlang::quos(...)
  dfnames <- sapply(quos_args, rlang::quo_name)

  dfnames1 <- names(list(...))
  if(!is.null(dfnames1)) dfnames[nchar(dfnames1)>0] <- dfnames1[nchar(dfnames1)>0]

  dflist <- lapply(quos_args, rlang::eval_tidy)
  dfs <- lapply(dflist, as.data.frame)

  for (i in which(sapply(dfs, ncol) == 1)){ colnames(dfs[[i]]) <- dfnames[i] }
  for (i in which(sapply(dfs, ncol) == 1)){ dfs[[i]] <- t(dfs[[i]]) }

  allids <- lapply(dfs, colnames)
  ids <- Reduce(f = intersect, x = allids)

  rows <- unlist(lapply(dfs, rownames))
  if (any(duplicated(rows))) rows <- paste0(rows, ".", rep(dfnames, sapply(dfs, nrow)))

  subdfs <- lapply(dfs, function(tmp) tmp[,ids,drop = FALSE] )
  res <- Reduce(f = rbind, x = subdfs)
  rownames(res) <- rows
  res
}














#' Write dataframes to .xlsx file
#'
#' @description
#' Write a dataframe or a named list of dataframes to a .xlsx file
#'
#' @param data
#' @param filename
#' @param adjwidths
#'
#' @export
#'
#' @seealso openxlsx::writeData
#' @examples
#' writeTables(mtcars, file = "example.xlsx")
#' writeTables(list(mtcars = mtcars, iris = iris), file = "example.xlsx")
#'
# writeTables <- function(data, file, rowNames = TRUE, adjwidths = TRUE, ...){
#
#   if (!"list" %in% class(data)){
#     data <- setNames(list(data), gsub("\\..*$", "", basename(file)))
#   }
#
#   newnames <- cutstr(names(data), maxchar = 29)
#   if (any(duplicated(newnames))){
#     newnames <- cutstr(names(data), maxchar = 26)
#     newnames <- dedupl(newnames)
#   }
#   names(data) <- newnames
#
#   wb <- openxlsx::createWorkbook()
#
#   invisible(lapply(names(data), function(tmpname){
#     tmpdata <- data[[tmpname]]
#     openxlsx::addWorksheet(wb, tmpname)
#     openxlsx::writeData(wb, sheet = tmpname, x = tmpdata, rowNames = rowNames, ...)
#     if (adjwidths == TRUE){
#       openxlsx::setColWidths(wb, sheet = tmpname, cols = 1:(ncol(tmpdata)+1), widths = "auto")
#     }
#   }))
#
#   if (baseext(file) != "xlsx") file <- paste0(file, ".xlsx")
#   openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)
# }






#' Write dataframes to .xlsx file
#'
#' @description
#' Write a dataframe or a named list of dataframes to a .xlsx file
#'
#' @param data dataframe or a named list of dataframes
#' @param file file
#' @param rowNames print rownames to sheet
#' @param adjwidths adjust column widths
#' @param header format headers in bold font
#' @param scale_styles color gradient styles, e.g. scale_styles = list(NES = c("#6200d0" = -1, "#ffffff" = 0, "#ffca0a" = 1))
#' @param condition_styles conditional highlighting styles, e.g. condition_styles = list(padj = c(fontColour = "#ff2222", rule = "<=0.05"))
#' @param highlight_styles value highlighting styles, e.g. highlight_styles = list(genes = list(fgFill = "#3be1ff", values = c("MTOR", "EGFR")))
#' @param check test if .xlsx file can be re-imported without changes to data
#' @param ...
#'
#' @return
#' @export
#' @seealso openxlsx::writeData
#' @examples
#' writeTables(mtcars, file = "example.xlsx")
#' writeTables(list(mtcars = mtcars, iris = iris), file = "example.xlsx")
#' writeTables(mtcars, file = "example.xlsx", scale_styles = list(mpg = c("#ffffff", "#6200d0"), cyl = c("#ffffff" = 0, "#ffca0a" = 10)))
#' writeTables(mtcars, file = "example.xlsx", condition_styles = list(mpg = c(fontColour = "#6200d0", rule = ">20")))
#' writeTables(iris, file = "example.xlsx", highlight_styles = list(Species = list(fgFill = "#3be1ff", values = c("setosa", "virginica"))), rowNames = FALSE)
writeTables <- function(data, file, rowNames = TRUE, adjwidths = TRUE,
                        header = "bold",
                        scale_styles = list(NES = c("#ffca0a" = 4, "#ffffff" = 0, "#6200d0" = -4)),
                        condition_styles = list(padj = c(fontColour = "#ff2222", bgFill = "#3be1ff", rule = "<=0.05")),
                        highlight_styles = list(term = list(fgFill = "#3be1ff", values = c())),
                        check = FALSE, ...){

  stopifnot(requireNamespace("openxlsx"))

  datamisc::colorcat("Add color gradients with 'scale_styles'", col = "blue")
  datamisc::colorcat("Add conditional highlighting with 'condition_styles'", col = "blue")
  datamisc::colorcat("Highlight certain values with 'highlight_styles'", col = "blue")
  datamisc::colorcat("Additional style arguments are passed to 'createStyle()'", col = "blue")


  # helper functions
  .get_scale_style <- function(style){
    if (is.null(names(style))){
      colors <- style
      limits <- NULL
    } else {
      limits <- sort(style)
      colors <- names(limits)

    }
    list(style = colors, rule = limits)
  }

  .get_condition_style <- function(style){
    style <- as.list(style)
    rule <- style$rule
    style$rule <- NULL
    style <- do.call(openxlsx::createStyle, style)
    list(style = style, rule = rule)
  }

  .get_highlight_style <- function(style){
    values <- style[["values"]]
    style["values"] <- NULL
    style <- do.call(openxlsx::createStyle, style)
    list(style = style, values = values)
  }


  if (!"list" %in% class(data)){
    data <- setNames(list(data), gsub("\\..*$", "", basename(file)))
  }

  # adjust sheet names
  newnames <- cutstr(names(data), maxchar = 29)
  if (any(duplicated(newnames))){
    newnames <- cutstr(names(data), maxchar = 26)
    newnames <- dedupl(newnames)
  }
  names(data) <- newnames


  # define styles
  if (!is.null(header)){
    header_style <- openxlsx::createStyle(textDecoration = header)
  }

  # create workbook
  wb <- openxlsx::createWorkbook()

  # add sheets
  table <- invisible(lapply(names(data), function(tmpname){

    tmpdata <- as.data.frame(data[[tmpname]])
    openxlsx::addWorksheet(wb, tmpname)

    # write data to sheet
    openxlsx::writeData(wb, sheet = tmpname, x = tmpdata, rowNames = rowNames, headerStyle = header_style, ...)

    if (adjwidths == TRUE){
      openxlsx::setColWidths(wb, sheet = tmpname, cols = 1:(ncol(tmpdata) + as.numeric(rowNames)), widths = "auto")
    }

    # apply styles: color scales
    if (length(scale_styles) > 0){
      for (col in names(scale_styles)){
        if (!col %in% colnames(tmpdata)) next
        col_style <- .get_scale_style(scale_styles[[col]])
        openxlsx::conditionalFormatting(wb,
                                        sheet = tmpname,
                                        cols = which(colnames(tmpdata) %in% col) + as.numeric(rowNames),
                                        rows = 1 + 1:nrow(tmpdata), # all rows
                                        rule = col_style$rule,
                                        style = col_style$style,
                                        type = "colourScale")
      }
    }

    # apply styles: conditional highlighting
    if (length(condition_styles) > 0){
      for (col in names(condition_styles)){
        if (!col %in% colnames(tmpdata)) next
        col_style <- .get_condition_style(condition_styles[[col]])
        openxlsx::conditionalFormatting(wb,
                                        sheet = tmpname,
                                        cols = which(colnames(tmpdata) %in% col) + as.numeric(rowNames),
                                        rows = 1 + 1:nrow(tmpdata), # all rows
                                        rule = col_style$rule,
                                        style = col_style$style)
      }
    }


    # apply styles: value highlighting
    if (length(highlight_styles) > 0){
      for (col in names(highlight_styles)){
        if (!col %in% colnames(tmpdata)) next
        col_style <- .get_highlight_style(highlight_styles[[col]])
        openxlsx::addStyle(wb,
                           sheet = tmpname,
                           cols = which(colnames(tmpdata) %in% col) + as.numeric(rowNames),
                           rows = 1 + which(tmpdata[[col]] %in% col_style$values), # select rows
                           style = col_style$style)
      }
    }

    tmpdata
  }))


  # save file
  if (baseext(file) != "xlsx") file <- paste0(file, ".xlsx")
  openxlsx::saveWorkbook(wb, file = file, overwrite = TRUE)

  # check if writing to excel changed any of the values
  if (check == TRUE){
    table_check <- readTables(file, rowNames = rowNames)
    if (is.data.frame(table_check)) table_check <- list(table_check)
    for (i in 1:length(table_check)){
      print(all.equal(table[[i]], table_check[[i]]))
    }
  }
}




#' Read all sheets from .xlsx file
#'
#' @param file
#' @param rowNames
#' @param ...
#'
#' @export
#'
#' @examples
#' writeTables(list(mtcars = mtcars, iris = iris), file = "example.xlsx")
#' readTables(file = "example.xlsx")
readTables <- function(file, rowNames = TRUE, ...){
  sheets <- openxlsx::getSheetNames(file)
  res <- lapply(setNames(sheets, sheets), function(tmp) openxlsx::read.xlsx(sheet = tmp, xlsxFile = file, rowNames = rowNames, ...) )

  # add a check for common excel gene fails


  if (length(res) == 1) res <- res[[1]]
  res
}








### Apply functions ------



#' Nested apply
#'
#' @param X List
#' @param FUN Function
#' @param n Level
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
napply <- function(X, FUN, n, ...){

  if (n == 0) return( FUN(X) )
  if (n == 1) return( lapply(X, FUN) )

  args <- list(...)
  args <- unlist(lapply(seq_along(args), function(i) paste0(names(args)[i], "=", as.character(args[[i]])) ))
  args <- paste(args, collapse = ", ")

  call <- paste("X", paste(rep("lapply", n-1), collapse = ", "), "FUN", args, sep = ", ")
  call <- paste0("lapply(", call, ")")

  eval(parse(text = call))
}







#' Drop-in replacement for lapply to quickly identify all failing items
#'
#' @param X
#' @param FUN
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
dbugapply <- function(X, FUN, ...){

  if (is.null(names(X)) | any(duplicated(names(X)))) stop("List must be uniquely named!")
  x <- names(X)

  res <- lapply(setNames(x, x), function(tmp){
    tryCatch(
      {
        FUN(X[[tmp]])
      },
      error = function(msg){
        message(paste0("Failed on ", tmp, ":"))
        message(msg)
      }
    )
  })

  invisible(res)
}





### Data processing ------



#' Hierarchical clustering of data
#'
#' @param data matrix or dataframe
#' @param method hclust
#' @param rows cluster rows
#' @param cols cluster columns
#' @param inf handling of non-finite values
#' @param na handling of missing values
#' @param ...
#'
#' @return clust
#' @export
#'
#' @examples
clusterData <- function(data, method = "hclust", rows = NULL, cols = NULL, inf = NULL, na = NULL, ...){

  ### Cluster rows/columns of a dataframe or matrix with NA or Inf values


  # input arguments
  if (is.null(rows)) rows <- nrow(data) < 1000
  if (is.null(cols)) cols <- ncol(data) < 1000
  if (is.null( na))  na <- ""
  if (is.null(inf)) inf <- ""



  # +/-Inf value handling
  if (any(!is.finite(nat(data)))){

    if (naf(inf == FALSE)) data <- subInf(data)
    if (is.na(inf)) data[!is.finite(data)] <- inf

  }



  # NA value handling (1)
  if (any(is.na(data))){


    if (is.numeric(na)) data[is.na(data)] <- na



  }







  # Clustering
  tmp <- list(rows = data, cols = t(data))
  res <- lapply(tmp[c(rows, cols)], function(data){

    # NA value handling (2)
    if (na == "omit") data <- na.omit(data)

    clust <- NULL

    if (method == "hclust") clust <- dendsort::dendsort(stats::hclust(stats::dist(data)))







    clust
  })



  # col.clusters <- dendextend::rotate(col.clusters, cluster_order)

  res
}




#' Make directory
#'
#' Recursively make a new directory (if not existing)
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
mkdir <- function(path){
  stopifnot(class(path) == "character")
  if (!dir.exists(path)){
    dir.create(path = path, recursive = TRUE)
  }
}





#' Scale data within range
#'
#' @param data Matrix or vector
#' @param from minimum
#' @param to maximum
#'
#' @return
#' @export
#'
#' @examples
#' rangescale(1:5)
rangescale <- function(data, from = 0, to = 1){

  data_min <- min(data, na.rm = TRUE)
  data_max <- max(data, na.rm = TRUE)

  (data - data_min)/(data_max - data_min) * (to - from) + from

}






#' Summarise multiple columns by grouping
#'
#' @param data
#' @param coldata
#' @param by
#' @param FUN
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summarise_cols <- function(data, coldata = NULL, by = NULL, FUN = NULL, ...){

  if (!is.null(coldata)){
    by <- rlang::enquo(by)
    coldata <- coldata[colnames(data),, drop = FALSE]
    grouping <- dplyr::pull(coldata, !!by)
  } else {
    grouping <- by
  }

  stopifnot(length(grouping) == ncol(data))

  groups <- unique(grouping)

  res <- as.data.frame(lapply(groups, function(g){
    FUN(data[, naf(grouping == g), drop = FALSE], ...)
  }))
  colnames(res) <- groups

  res
}


#' Summarise multiple rows by grouping
#'
#' @param data
#' @param rowdata
#' @param by
#' @param FUN
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
summarise_rows <- function(data, rowdata = NULL, by = NULL, FUN = NULL, ...){

  if (!is.null(rowdata)){
    by <- rlang::enquo(by)
    rowdata <- rowdata[rownames(data),, drop = FALSE]
    grouping <- dplyr::pull(rowdata, !!by)
  } else {
    grouping <- by
  }

  stopifnot(length(grouping) == nrow(data))

  groups <- unique(grouping)

  res <- as.data.frame(lapply(groups, function(g){
    FUN(data[naf(grouping == g),, drop = FALSE], ...)
  }))
  rownames(res) <- groups

  res
}



#' Set the Names in an Object
#'
#' @param object
#' @param colnames
#'
#' @return
#' @export
#'
#' @examples
setColnames <- function(object, colnames = NULL){
  colnames(object) <- colnames
  object
}


#' Set the Names in an Object
#'
#' @param object
#' @param rownames
#'
#' @return
#' @export
#'
#' @examples
setRownames <- function(object, rownames = NULL){
  rownames(object) <- rownames
  object
}





#' Source R code of an RMD file
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
sourceRMD <- function(file){
  tmp <- tempfile()
  on.exit(unlink(tmp))
  invisible(knitr::knit(input = file, quiet = TRUE, output = tmp))
}








#' Format p-values for printing
#'
#' @param p
#' @param min
#' @param scientific
#' @param add
#' @param stars
#'
#' @return
#' @export
#'
#' @examples
pval_format <- function(p, min = 0.001, scientific = NULL, add = "p ", stars = FALSE){

  if (!is.null(min)){
    digits <- nchar(sub("^-?\\d*\\.?","", min))
    ptext <- round(p, digits = digits)
    poob <- p < min
    ptext[poob] <- min
    ptext <- paste0(add, ifelse(poob, "<", "="), " ", ptext)

    psc <- signif(p, digits) |> format(scientific=TRUE)
    if (!is.null(scientific)) ptext[poob] <- psc[poob]

  }

  # format(p, scientific=TRUE)

  ptext
}







#' Combine character values into a vector or list without quotes
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' ce(a,b,c)
ce <- function(...){
  args <- rlang::enquos(...)
  sapply(args, rlang::as_name)
}



#' Source R functions from dir into local env
#'
#' @param dir
#' @param pattern
#'
#' @return
#' @export
#'
#' @examples
#' newenv <- sourcelib("lib")
#' attach(newenv)
sourcelib <- function(dir = "lib", pattern = "\\.r$|\\.R$"){
  env <- new.env()
  grep(pattern = ".r$", list.files(dir, full.names = TRUE), value = TRUE) |> sapply(source, local = env) |> invisible()
  env
}




#' Get data from the GEO database
#'
#' @param ID
#' @param geodir
#' @param fetch_files
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' GSE234099 <- getGEOdata(ID = "GSE234099", fetch_files = FALSE)
#' GSE222693 <- getGEOdata(ID = "GSE222693")
getGEOdata <- function(ID = "GSE222693", geodir = "~/myScratch/GEO", fetch_files = TRUE, ...){

  cat(crayon::blue("Use 'Biobase::pData($geo$series_matrix.txt.gz)' to get the samplesheet\n"))
  cat(crayon::blue("Use 'Biobase::exprs($geo$series_matrix.txt.gz)' to get the expression data (if present)\n"))

  stopifnot(requireNamespace("GEOquery"))

  dir.create(geodir, showWarnings = FALSE)
  subdir <- file.path(geodir, ID)
  dir.create(subdir, showWarnings = FALSE)
  if (length(list.files(subdir)) > 0){
    warning("Warning: Directory not empty!")
  }

  geo <- GEOquery::getGEO(ID, destdir = subdir, ...)
  supp <- GEOquery::getGEOSuppFiles(makeDirectory = FALSE, "GSE222693", baseDir = subdir, fetch_files = fetch_files)

  if (fetch_files == FALSE){
    suppfiles <- unname(supp["url"])
  } else {
    suppfiles <- rownames(supp)
  }

  list(geo = geo, suppfiles = suppfiles)
}



#' Safe filenames
#'
#' @param names
#' @param replacement
#' @param unique
#' @param sep
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' make_filenames("file/na.me")
make_filenames <- function(names, replacement = "_", unique = TRUE, sep = .Platform$file.sep, ...){
  names <- sub(sep, "_", names, fixed = TRUE)
  names <- make.names(names, unique = unique, ...)
  names <- sub(".", replacement, names, fixed = TRUE)
  names <- sub(".", replacement, names, fixed = TRUE)
  names
}




#' Software version reporting
#'
#' @return
#' @export
#'
#' @examples
versionInfo <- function(){

  env <- sessioninfo::platform_info()
  df1 <- t(as.data.frame(env)[,c("version", "os")])
  df1 <- cbind(c("R", "OS"), df1)
  df1[1,2] <- gsub("R version ", "", df1[1,2])

  pkgs <- sessioninfo::package_info(pkgs = .packages(), dependencies = FALSE)
  df2 <- as.data.frame(pkgs)[,c("package", "loadedversion", "source")]
  colnames(df2) <- c("package", "version", "source")

  list(platform = df1, packages = df2)
}




#' Intersection of all elements
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' intersection(LETTERS[1:10], LETTERS[5:8], LETTERS[7:12])
intersection <- function(...){

  args <- list(...)
  if (length(args) == 1 & is.list(args[[1]])){
    args <- args[[1]]
  }

  Reduce(f = intersect, x = args)

}



















