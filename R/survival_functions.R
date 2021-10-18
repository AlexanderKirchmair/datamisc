

#' Stratify variables based on quantiles
#'
#' @param x Data vector
#' @param n Number of strata
#' @param ties Assign quantile ties to higher or lower group
#'
#' @return
#' @export
#'
#' @examples
stratify <- function(x, n = 2, ties = "lower", optimal = FALSE, surv = NULL, ...){

  if (optimal == FALSE){

    if ("data.frame" %in% class(x)){
      ix <- sapply(x, is.numeric)
      ix[colnames(x) %in% c("surv")] <- FALSE

      x_strat <- as.data.frame(lapply(x[,ix,drop=FALSE], function(y){
        stratify_var(x = y, n = n, ties = ties)
      }))
      x[,ix] <- x_strat

    } else {
      x <- stratify_var(x = x, n = n, ties = ties, ...)
    }

  } else {

    if (is.null(dim(x))) x <- data.frame(x = x)
    if (is.null(surv)) surv <- x$surv
    stopifnot(!is.null(surv))

    ix <- sapply(x, is.numeric)
    ix[colnames(x) %in% c("surv")] <- FALSE

    tmp <- data.frame(x[,ix,drop = FALSE],
                      time = surv[,1],
                      status = surv[,2])

    cutpoints <- survminer::surv_cutpoint(tmp, time = "time", event = "status", variables = setdiff(colnames(tmp), c("time", "status")),
                                          progressbar = FALSE, ...)
    cutdata <- select(data.frame(survminer::surv_categorize(x = cutpoints)), -c(time, status))

    tmp <- as.data.frame(lapply(setNames(colnames(cutdata), colnames(cutdata)), function(var){
      val <- signif(cutpoints$cutpoint[var, "cutpoint"], 2)
      stratvar <- ifelse(cutdata[[var]] == "low", paste0("<", val), paste0(">", val))
      relevel(factor(stratvar), ref = paste0("<", val))
    }))

    x[,ix] <- tmp
    if (ncol(x) == 1) x <- x[[1]]
  }

  x
}



#' Stratify helper function
#'
stratify_var <- function(x, n = 2, ties = "lower", ...){

  if (n == 0) return(rep(mean(x), n = length(x)))
  if (n == 1) return(x)

  # get quantiles
  qs <- quantile(x, seq_along(1:(n-1))/n, na.rm = TRUE)

  # get quantile names
  if (n == 2){
    classes.names <- c("Lo", "Hi")
  } else if (n == 3) {
    classes.names <- c("Lo", "Mid", "Hi")
  } else {
    names(qs) <- paste0( as.character(round( 100*seq_along(1:(n-1))/n, 1)), "%")
    tmpnames <- c()
    if (n > 2){
      for (i in 1:(length(qs)-1)){
        tmpnames[i] <- paste0(names(qs)[i], "-", names(qs)[i+1])
      }
    }
    classes.names <- c(paste0("<", names(qs)[1]), tmpnames, paste0(">", names(qs)[length(qs)]))
  }

  classes <- rep(classes.names[1], length(x))

  if (ties == "higher"){
    for (i in 1:length(qs)){
      classes[ x >= qs[i] ] <- classes.names[i+1]
    }
  } else if (ties == "lower"){
    for (i in 1:length(qs)){
      classes[ x > qs[i] ] <- classes.names[i+1]
    }
  }

  if (n == 2){
    classes <- relevel(factor(classes, levels = classes.names), ref = "Lo")
  } else {
    classes <- factor(classes, levels = classes.names, ordered = TRUE)
  }

  names(classes) <- names(x)
  classes[is.na(x)] <- NA

  return(classes)
}


