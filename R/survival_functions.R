

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



















ggplotForest <- function(df, x = H, y = term, color = -log10(pval), ptres = 0.05, ...){

  # Input
  df <- as.data.frame(df)

  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  color <- rlang::enquo(color)
  xcol <- rlang::as_name(x)
  ycol <- rlang::as_name(y)

  siglab <- c(sig = "pval < 0.05", nonsig = "pval > 0.05")

  df <- dplyr::arrange(df, pval)

  if (is.null(df[[ycol]])){
    df[[ycol]] <- factor(rownames(df), ordered = TRUE, levels = rev(rownames(df)))
  } else {
    df[[ycol]] <- factor(df[[ycol]], ordered = TRUE, levels = rev(as.character(df[[ycol]])))
  }
  df$sig <- relevel(factor(ifelse(df$pval <= ptres, siglab["sig"], siglab["nonsig"])), ref = siglab["nonsig"])
  df$odd <- ifelse(1:nrow(df) %% 2 == 1, "odd", "even")

  grid_line <- element_line(colour = "grey70", size = 0.3, linetype = "dashed")


  # Adjust limits
  xmin <- min(df$CI95_lower, na.rm = TRUE)
  xmin <- floor(xmin * 10^-round(log10(xmin))) / 10^-round(log10(xmin)) # round down
  xmax <- max(df$CI95_upper, na.rm = TRUE)
  xmax <- ceiling(xmax * 10^-round(log10(xmax))) / 10^-round(log10(xmax)) # round up
  xlim <- c(xmin, xmax)
  df$CI95_lower[df$CI95_lower < xlim[1]] <- xlim[1]
  df$CI95_upper[df$CI95_upper > xlim[2]] <- xlim[2]
  # add v. segments if cut-off

  # Ggplot
  ggplot(df, aes(x = !!x, y = !!y, color = !!color)) +
    theme_minimal() +
    theme(axis.text = element_text(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.x = element_line(colour = "black"),
          panel.grid.major.x = grid_line,
          panel.grid.major.y = element_blank(),
          panel.grid.minor = grid_line) +
    geom_rect(aes(fill = odd), color = NA, show.legend = FALSE,
              ymin = as.numeric(df[[ycol]])-0.5, ymax = as.numeric(df[[ycol]])+0.5, xmin = log10(xlim[1]), xmax = log10(xlim[2])) +
    geom_vline(xintercept = 1, colour = "grey50") +
    geom_segment(aes(x = CI95_lower, xend = CI95_upper, y = !!y, yend = !!y), color = "black") +
    geom_point(aes(shape = sig), size = 2.8) +
    scale_x_continuous(limits = xlim, trans = "log10", expand = ggplot2::expansion(0,0)) +
    xlab("") +
    coord_cartesian(clip = "off") +
    scale_fill_manual(values = c("odd" = rgb(1,1,1,0), "even" = rgb(0,0,0,0.1))) +
    scale_color_gradient(low = rgb(0,0,0), high = rgb(1,0,0)) +
    scale_shape_manual("signif.", values = setNames(c(0,15), siglab[c("nonsig", "sig")])) +
    ylab("")

}




results <- function(x, ...) {
  UseMethod("results", x)
}

results.survfit <- function(x, ...){

  km <- x
  data <- km$data
  survdiff.call <- str2lang(paste0("survdiff(", deparse(km$formula), ", data = data)"))
  diff <- eval(survdiff.call)
  diff$pval <- 1 - pchisq(diff$chisq, length(diff$n)-1)

  diff
}




KM <- function(...){

  args <- list(...)
  args.class <- lapply(args, class)

  ix <- c("formula" = which(sapply(args.class, function(x) "formula" %in% x)),
          "data" = which(sapply(args.class, function(x) "data.frame" %in% x | "data.matrix" %in% x | "matrix" %in% x)),
          "group" = which(sapply(args.class, function(x) "factor" %in% x | "character" %in% x)))

  if ("formula" %in% names(ix)) formula <- args[[ ix["formula"] ]] else formula <- NULL
  if ("data" %in% names(ix)) data <- args[[ ix["data"] ]]  else data <- NULL
  if ("group" %in% names(ix)) group <- args[[ ix["group"] ]]  else group <- NULL

  if (is.null(formula)){
    data$surv <- survival::Surv(time = data[,"time"], event = data[,"status"])

    if (!is.null(group)){
      stopifnot(length(group) == nrow(data))
      formula <- surv ~ group
    } else {
      formula <- as.formula(paste0("surv ~ ", paste(setdiff(colnames(data), c("time", "status", "surv")), collapse = " + ")))
    }

  }

  survfit.call <- str2lang(paste0("survfit(", deparse(formula), ", data = data)"))
  km <- eval(survfit.call)
  km$data <- data
  km$formula <- formula

  if (length(labels(terms(as.formula(km$call$formula)))) > 0){
    survdiff.call <- str2lang(paste0("survdiff(", deparse(km$formula), ", data = data)"))
    km$diff <- eval(survdiff.call)
    km$pval <- 1 - pchisq(km$diff$chisq, length(km$diff$n)-1)
  }

  km
}









ggplotSurv <- function(fit, data = NULL, xlabel = "time", ylabel = "survival", title = NULL, class_title = "", rename = TRUE, censorsize = 3, colors = NULL, CI = FALSE, base_size = 18, ...){

  # Function to plot survival curves

  if (!is.null(data)) fit$data <- recodeCoxdata(data)
  stopifnot(!is.null(fit$data))

  if ("coxph" %in% class(fit)){

    coxfit <- fit
    vars <- labels(terms(as.formula(coxfit$call$formula)))
    tmp <- vars[grep("\\(*\\)", vars)]
    tmp <- gsub(".*\\(|\\).*", "", tmp)
    vars[grep("\\(*\\)", vars)] <- tmp


    newdf <- coxfit$data[,vars,drop = FALSE]
    newdf <- as.data.frame(lapply(newdf, function(x) if (is.numeric(x)) stratify(x) else x))
    newdf <- unique(newdf)
    newdf <- newdf[rowSums(is.na(newdf)) == 0,,drop=FALSE]

    newdf2 <- as.data.frame(lapply(colnames(newdf), function(col) paste0(col, "-", as.character(newdf[[col]])) ))
    newgroups <- apply(newdf2, 1, function(row) paste0(as.character(row), collapse = "/"))
    newgroups <- gsub("__", "", newgroups)
    rownames(newdf) <- newgroups

    fit <- survfit(coxfit, newdata = newdf)
    fit$data <- coxfit$data
    fit$groups <- newgroups # 'fit$strata' not compatible with ggsurvplot

  }


  if (!is.null(fit$strata)) fit$groups <- names(fit$strata)
  vars <- labels(terms(fit$formula))
  if (is.null(title)) title <- paste(vars, collapse = " - ")

  pchar <- function(p, digits = 3, min_p = 0.001, str = "p"){
    if (p < min_p){
      str <- paste0(str, " < ", min_p)
    } else {
      str <- paste0(str, " = ", round(p, digits))
    }
    str
  }

  if (is.null(colors)){
    if (length(fit$groups) <= 8){
      colors <- as.character(paletteer::paletteer_d("colorblindr::OkabeIto")[c(5,1,3,7,4,2,8,6)][seq(fit$groups)])
    } else if (length(fit$groups) <= 256){
      colors <- as.character(paletteer::paletteer_d("palettesForR::Coldfire")[floor(seq(1, 256, length = length(fit$groups)))])
    } else {
      colors <- circlize::rand_color(length(fit$groups))
    }
  } else {
    colors <- colors[seq(fit$groups)]
  }

  if (length(vars) == 1){

    if (!is.null(names(colors))){
      gnames <- names(colors)
      names(colors) <- paste0(vars, "=", names(colors))
    } else {
      if ("factor" %in% class(fit$data[[vars]])){
        gnames <- levels(fit$data[[vars]])
        names(colors) <- paste0(vars, "=", gnames)
      } else {
        gnames <- fit$groups # sort(unique(fit$data[[vars]]))
        names(colors) <- gnames
        gnames <- sub(".*=", "", gnames)
      }

    }

  }

  pstr <- NULL
  if (!is.null(fit$pval)) pstr <- pchar(p = fit$pval)
  if (!is.null(fit$padj)) pstr <- paste0(pstr, "\n", pchar(p = fit$padj, str = "adj"), "")

  gp <- ggsurvplot(fit = fit,
                   data = fit$data,
                   risk.table = FALSE,
                   pval = pstr,
                   # palette = setNames(colors, NA),
                   ggtheme = theme_classic(base_size = base_size),
                   pval.size = base_size*0.35,
                   risk.table.y.text.col = TRUE,
                   risk.table.y.text = FALSE,
                   xlab = xlabel,
                   ylab = ylabel,
                   censor.shape = 124,
                   censor.size = censorsize,
                   conf.int = ifelse(CI, TRUE, FALSE),
                   legend.title = class_title)

  suppressMessages(
    gp$plot <- gp$plot +
      scale_x_continuous(expand = c(0, 0.2)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
      theme(axis.text = element_text(colour = "black", size = base_size),
            axis.line = element_line(lineend = "square"),
            legend.position = "right",
            legend.background = element_blank(),
            legend.text = element_text(size = rel(0.7)),
            plot.margin = margin(0.05,0.1,0.05,0.05, unit = "npc"),
            plot.title = element_text(hjust = 0.5, size = base_size)) +
      coord_cartesian(clip = "off") +
      ggtitle(title))

  if (rename == TRUE) gp$plot <- gp$plot + scale_color_manual(values = colors, labels = gnames)

  return(gp$plot)
}





quo_class <- function(x){
  x_class <- NA
  try({x_class <- class(rlang::eval_tidy(x))}, silent = TRUE)
  x_class
}




COX <- function(formula, data = NULL, update = NULL, surv = NULL, robust = FALSE, response = "surv", ...){

  if ("character" %in% class(formula)){
    formula <- reformulate(formula, response = response)
  }

  if ("coxph" %in% class(formula)){
    if (is.null(data)) data <- formula$data
    formula <- formula$formula
  }

  formula2 <- rlang::enquo(update)

  if (!rlang::quo_is_null(formula2)){
    if (all(quo_class(formula2) == "formula")){
      formula2 <- update
    } else {
      fstr <- rlang::as_name(formula2)

      formula2 <- as.formula(paste0("~ . + ", paste(fstr, collapse = " + ")))
    }

    formula <- update.formula(formula, formula2)
  }

  data <- recodeCoxdata(data)

  if (!robust) cox.call <- str2lang(paste0("survival::coxph(", deparse1(formula), ", data = data, ...)"))
  if ( robust) cox.call <- str2lang(paste0("coxrobust::coxr(", deparse1(formula), ", data = data, ...)"))

  coxfit <- eval(cox.call)
  coxfit$data <- data

  if (!is.null(rownames(data))) coxfit$check <- eval(str2lang(paste0("survival::survcheck(", deparse1(formula), ", data = data, id = rownames(data), ...)")))

  return(coxfit)
}




results <- function(x, ...) {
  UseMethod("results", x)
}

results.survfit <- function(x, ...){

  km <- x
  data <- km$data
  survdiff.call <- str2lang(paste0("survdiff(", deparse(km$formula), ", data = data)"))
  diff <- eval(survdiff.call)
  diff$pval <- 1 - pchisq(diff$chisq, length(diff$n)-1)

  diff
}


results.coxph <- function(x, zph = TRUE, vif = TRUE, sort = TRUE, ...){

  tmp <- summary(x)

  if ("coxph.penal" %in% class(x)){

    varl <- strsplit(rownames(tmp$coefficients), split = ":")
    res.df <- data.frame(var = sapply(varl, function(v) paste(sub("__.*", "", v), collapse = ":")),
                         level = sapply(varl, function(v) paste(sub(".*__", "", v), collapse = ":")),
                         H = exp(tmp$coefficients[,"coef"]),
                         pval = tmp$coefficients[,"p"],
                         CI95_lower = data.frame(tmp$conf.int, check.names = F)[rownames(tmp$coefficients),][,"lower .95"],
                         CI95_upper = data.frame(tmp$conf.int, check.names = F)[rownames(tmp$coefficients),][,"upper .95"])

    res.df$n <- tmp$n
    res.df$LRT.model <- tmp$logtest["pvalue"]

    if (!any(duplicated(rownames(tmp$coefficients)))) rownames(res.df) <- rownames(tmp$coefficients)


  } else {

    varl <- strsplit(rownames(tmp$coefficients), split = ":")
    res.df <- data.frame(var = sapply(varl, function(v) paste(sub("__.*", "", v), collapse = ":")),
                         level = sapply(varl, function(v) paste(sub(".*__", "", v), collapse = ":")),
                         H = tmp$coefficients[,"exp(coef)"],
                         pval = tmp$coefficients[,"Pr(>|z|)"],
                         CI95_lower = tmp$conf.int[,"lower .95"],
                         CI95_upper = tmp$conf.int[,"upper .95"],
                         n = tmp$n,
                         LRT.model = tmp$logtest["pvalue"],
                         Wald.model = tmp$waldtest["pvalue"],
                         Sctest.model = tmp$sctest["pvalue"])

    rownames(res.df) <- rownames(tmp$coefficients)
  }


  res.df$AIC.model <- extractAIC(x)[2]
  res.df$VIF <- rms::vif(x)[rownames(res.df)]
  try({if (zph == TRUE) res.df$zph <-  data.frame(survival::cox.zph(x)$table)[res.df$var,"p"]})

  if (sort == TRUE) res.df <- arrange(res.df, pval)

  return(res.df)
}







results.list <- function(x, ...){

  lapply(x, results, ...)

}



coxtable <- function(coxlist, vars = NULL, padj = "holm", filter = TRUE, ...){

  if (!is.null(names(coxlist))){
    names(coxlist) <- gsub(" * ", ":", names(coxlist), fixed = TRUE)
    names(coxlist) <- gsub("* ", ":", names(coxlist), fixed = TRUE)
    names(coxlist) <- gsub(" *", ":", names(coxlist), fixed = TRUE)
    names(coxlist) <- gsub(" : ", ":", names(coxlist), fixed = TRUE)
    names(coxlist) <- gsub(": ", ":", names(coxlist), fixed = TRUE)
    names(coxlist) <- gsub(" :", ":", names(coxlist), fixed = TRUE)
  }

  reslist <- results(coxlist, ...)
  if (is.null(names(reslist))) names(reslist) <- seq(reslist)
  reslist <- lapply(setNames(seq(reslist), names(reslist)), function(i){
    reslist[[i]]$model <- trimws(gsub(".*~", "", deparse1(formula(coxlist[[i]]$terms)))); reslist[[i]] })
  res <- Reduce(f = rbind, reslist)
  rownames(res) <- NULL

  if (filter == TRUE){
    if (is.null(vars)) vars <- unique(unlist(strsplit(names(reslist), " ")))

    varl2 <- strsplit(unique(res$var), ":")
    vars <- sapply(strsplit(vars, ":"), function(x){
      unique(res$var)[ sapply(varl2, function(v2) setequal(x, v2)) ]
    })

    # vars <- vars[vars %in% unique(res$var)]

    resdf <- res[res$var %in% vars,,drop = FALSE]
    rownames(resdf) <- NULL
  } else {
    resdf <- res
  }


  resdf$padj <- p.adjust(resdf$pval, method = padj)
  resdf <- dplyr::relocate(.data = resdf, "padj", .after = "pval")

  return(resdf)
}






selectCox <- function(full, reduced = NULL, FUN = AIC, ptres = 0.05, dHtres = 0.1){

  # full: all variables
  # reduced: variables to keep (subset of full)

  cox0 <- full
  val0 <- FUN(cox0)
  vars <- labels(terms(cox0))
  if (!is.null(reduced)) keep <- labels(terms(reduced)) else keep <- c()
  keep <- union(keep, grep("\\(.*\\)", vars, value = TRUE))
  if (length(setdiff(vars, keep)) == 0) return(full)

  res <- results(cox0, zph = FALSE)
  p <- setNames(res$pval, res$var)
  nonsigvars <- names(p[naf(p > ptres)])
  nonsigvars <- nonsigvars[!duplicated(nonsigvars)]
  sigvars <- names(p[naf(p <= ptres)])
  sigvars <- sigvars[!duplicated(sigvars)]

  # select vars to remove
  fits1 <- setdiff(vars, keep) %L>% function(x) COX(cox0, update = as.formula(paste0("~ . -", x)))
  val1 <- (fits1 %S>% function(x){FUN(x)}) %>% sort()
  remove <- names(val1[val1 < val0])
  remove <- setdiff(c(remove, setdiff(nonsigvars, remove)), keep)

  vifs <- rms::vif(cox0)
  vifs <- vifs[names(vifs) %in% remove] %>% sort(decreasing = TRUE) %>% names()

  # 1. Remove from full model
  cox1 <- list()
  cox1$a <- coxselectFUN(cox0, c(vifs, setdiff(remove, vifs)), FUN = FUN, rm = TRUE)
  cox1$b <- coxselectFUN(cox0, remove, FUN = FUN, rm = TRUE)
  cox1$c <- coxselectFUN(cox0, nonsigvars, FUN = FUN, rm = TRUE)
  cox1$d <- coxselectFUN(cox0, setdiff(vars, keep), FUN = FUN, rm = TRUE)

  cox1vars <- sort(table(unlist(lapply(cox1, function(tmp) labels(terms(tmp))))))
  d <- c(setdiff(remove, names(cox1vars)), setdiff(vars, names(cox1vars)))
  cox1 <- coxselectFUN(cox0, setdiff(d, keep), FUN = FUN, rm = TRUE)


  # 2. Re-add removed variables that alter estimates for the included variables
  removed <- setdiff(vars, labels(terms(cox1)))
  res <- results(cox1, zph = FALSE)
  H <- setNames(res$H, rownames(res))
  Hres <- removed %S>% function(var){
    res <- results(COX(cox1, update = as.formula(paste0("~ . +", var))), zph = FALSE, vif = FALSE)
    setNames(res$H, rownames(res))[names(H)]
  }
  mat <- abs(1 - Hres/H)
  mat[is.na(mat)] <- 0
  readd <- sort(setNames(colMaxs(mat), colnames(Hres)), decreasing = TRUE)
  readd <- readd[readd > dHtres]

  cox2 <- coxselectFUN(cox1, names(readd), FUN = FUN, rm = FALSE)


  # 3. Multiple selection
  cox <- cox2
  coxold <- ~ 1
  fin <- FALSE
  i <- 0
  while (!fin){

    if (setequal(labels(terms(cox)), labels(terms(coxold)))){
      if (FUN(cox) == FUN(coxold)){
        fin <- TRUE }}

    coxold <- cox
    i <- i + 1
    print(i)

    # drop all nonsig
    sigvars <- unique(subset(results(cox, zph = FALSE), pval <= ptres)$var)
    rmvars <- setdiff(labels(terms(cox)), c(sigvars, keep))
    cox <- COX(cox, update = as.formula(paste0("~ . -", paste0(rmvars, collapse = " - "))))

    # add1
    add1vars <- setdiff(rmvars, labels(terms(cox)))
    fits1 <- add1vars %L>% function(x) COX(cox, update = as.formula(paste0("~ . + ", x)))
    allsig <- sapply(fits1, FUN)

    # readd
    readd <- names(allsig[allsig < FUN(cox)])
    if (length(readd) > 0) cox<- COX(cox, update =  reformulate(c(".", readd)))
  }

  cox

}




elasticCOX <- function(cox, ...){

  vars <- labels(terms(cox$formula))
  stratvars <- gsub(".*\\(|)$", "", vars[grep("strata", vars)])

  svar <- setdiff(all.vars(cox$formula), c(vars, stratvars))
  vars <- vars[!grepl("strata", vars)]

  x <- cox$data[,vars]
  y <- cox$data[,svar]

  ix <- as.numeric(y[,1]) > 0 & rowSums(is.na(x)) == 0

  x <- subset(x, ix)
  y <- subset(y, ix)

  if (length(stratvars) > 0) for (sv in stratvars){ y <- stratifySurv(y, strata = cox$data[ix,sv]) }

  res <- glmnet(x, y, family = "cox", ...)

  # cv.glmnet(x, y, family = "cox") # only numeric?
  # coef(cvfit, s = "lambda.min")

  res
}

vars.coxph <- function(x, ...){
  labels(terms(x))
}


coxselectFUN <- function(cox, vars, FUN, rm = FALSE){
  if (rm == TRUE) con <- "-" else con <- "+"
  cox_new <- cox
  for (var in vars){
    cox_tmp <- COX(cox_new, update = as.formula(paste0("~ . ", con, var)))
    if (naf(FUN(cox_tmp) <= FUN(cox))){ cox_new <- cox_tmp }
  }
  cox_new
}






recodeCoxdata <- function(data, pattern = "__"){

  rows <- rownames(data)
  data <- as.data.frame(data)

  data <- lapply(data, function(tmp){

    if (is.factor(tmp)){
      tmp <- plyr::mapvalues(tmp, from = levels(tmp), to = sub(paste0("^", pattern), "", levels(tmp)))
      tmp <- plyr::mapvalues(tmp, from = levels(tmp), to = paste0("__", levels(tmp)))

    } else if (is.character(tmp)){
      tmp <- plyr::mapvalues(tmp, from = unique(tmp), to = sub(paste0("^", pattern), "", unique(tmp)))
      tmp <- plyr::mapvalues(tmp, from = unique(tmp), to = paste0("__", unique(tmp)))

    }
    tmp
  })

  data <- data.frame(row.names = rows, data)

  data
}


