



#' List pipe
#'
#' @param lhs List or list-like object
#' @param rhs Function
#'
#' @export
#'
#' @examples
`%L>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- listpipe(lhs = lhs, rhs = rhs, env = pe)
  return(res)
}



#' List assignment pipe
#'
#' @inheritParams `%L>%`
#' @export
#'
`%<L>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- invisible(listpipe(lhs = lhs, rhs = rhs, env = pe))
  assign(deparse(lhs), res, envir = pe)
}

#' List pipe
#'
#' @inheritParams `%L>%`
#' @export
#'
`%S>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- listpipe(lhs = lhs, rhs = rhs, env = pe, simplify = TRUE)
  return(res)
}

#' List assignment pipe
#'
#' @inheritParams `%S>%`
#' @export
#'
`%<S>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- invisible(listpipe(lhs = lhs, rhs = rhs, env = pe, simplify = TRUE))
  assign(deparse(lhs), res, envir = pe)
}

#' Parallel list pipe
#'
#' @inheritParams `%L>%`
#' @export
#'
`%P>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- listpipe(lhs = lhs, rhs = rhs, env = pe, parallel = TRUE)
  return(res)
}

#' Parallel list assignment pipe
#'
#' @inheritParams `%P>%`
#' @export
#'
`%<P>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- invisible(listpipe(lhs = lhs, rhs = rhs, env = pe, parallel = TRUE))
  assign(deparse(lhs), res, envir = pe)
}

#' SGE qsub list pipe
#'
#' @inheritParams `%L>%`
#' @export
#'
`%Q>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- listpipe(lhs = lhs, rhs = rhs, env = pe, qsub = TRUE)
  return(res)
}

#' SGE qsub list assignment pipe
#'
#' @inheritParams `%Q>%`
#' @export
#'
`%<Q>%` <- function(lhs, rhs){
  lhs <- substitute(lhs)
  rhs <- substitute(rhs)
  pe <- parent.frame()
  res <- invisible(listpipe(lhs = lhs, rhs = rhs, env = pe, qsub = TRUE))
  assign(deparse(lhs), res, envir = pe)
}



#' List pipe function
#'
#' @param lhs
#' @param rhs
#' @param env
#' @param parallel
#' @param simplify
#' @param qsub
#' @param ...
#'
#'
listpipe <- function(lhs, rhs, env, parallel = FALSE, simplify = FALSE, qsub = FALSE, ...){

  ### LHS (list) ----

  lhs <- eval(lhs, envir = env)
  if (is.null(names(lhs)) & (is.character(lhs) | is.integer(lhs))) names(lhs) <- as.character(lhs)

  ### RHS (function) ----

  FUN <- rlang::call_standardise(rhs)
  args <- rlang::call_args(FUN)

  # args are different for an anonymous function
  anon_f <- any(unlist(sapply(args, grepl, pattern = "function")))
  if (anon_f){
    FUN <- args[[length(args)]]
    FUN <- as.call(parse(text = as.character(FUN)))
    args <- unlist(args[-(length(args)-0:1)], recursive = FALSE)
  }

  # check if dot is specified in function args
  if (length(args) > 0){

    names(args)[names(args)==""] <- NA
    dix <- which(sapply(args, function(tmp) all(tmp == ".") ))

    # if dot is named, use name(s) later
    dot_args <- names(args)[dix]
    dot_args <- dot_args[!is.na(dot_args)]
    if (length(dot_args) == 0) dot_args <- NULL

    # if dot is unnamed, remove
    drm <- dix[is.na(names(args)[dix]) | is.null(names(args)[dix])]
    if (any(drm)) args <- args[-drm]
    if (length(args) == 0) args <- NULL

  } else dot_args <- NULL

  # remove args from function call (unless they are unnamed and passed to ...)
  args <- args[!is.na(names(args))]
  rm_args <- lapply(setNames(names(args), names(args)), function(tmp) rlang::zap() )
  FUN <- rlang::call_modify(FUN, ... = rlang::zap())
  FUN <- rlang::call_modify(FUN, !!!rm_args)

  # remove dots from call
  ftxt <- as.character(FUN)
  ftxt <- ftxt[ftxt != "."]
  FUN <- as.call(parse(text = ftxt))

  # remove dot from args
  if (length(dot_args) != 0) args[dot_args] <- NULL


  ### Iterate function over list ----

  itfun <- function(tmp, dot_args, args, FUN, env){
    # add tmp as dot or unnamed arg
    tmp_args <- list(tmp)
    if (!is.null(dot_args)) tmp_args <- rep(tmp_args, length(dot_args))
    all_args <- c(setNames(tmp_args, dot_args), args) # update current arguments
    FUN <- rlang::call_modify(FUN, !!!all_args) # add current arguments
    eval(FUN, envir = env)
  }


  if (parallel == TRUE){
    # parallel setup...
    if(BiocParallel::bpparam()$workers > 10) BiocParallel::register(BiocParallel::MulticoreParam(workers = 10))
    res <- BiocParallel::bplapply(X = setNames(seq_along(lhs), names(lhs)), FUN = function(i){
      assign(x = ".name", value = names(lhs)[i], envir = env)
      itfun(lhs[[i]], dot_args, args, FUN, env)
    }, BPPARAM = BiocParallel::bpparam(), ...)

  } else if (simplify == TRUE){
    res <- sapply(setNames(seq_along(lhs), names(lhs)), function(i){
      assign(x = ".name", value = names(lhs)[i], envir = env)
      itfun(lhs[[i]], dot_args, args, FUN, env)
    }, ...)

  } else if (qsub == TRUE){
    res <- qapply(setNames(seq_along(lhs), names(lhs)), env = env, function(i){
      assign(x = ".name", value = names(lhs)[i], envir = env)
      itfun(lhs[[i]], dot_args, args, FUN, env)
    })

  } else {
    res <- lapply(setNames(seq_along(lhs), names(lhs)), function(i){
      assign(x = ".name", value = names(lhs)[i], envir = env)
      itfun(lhs[[i]], dot_args, args, FUN, env)
    }, ...)

  }

  suppressWarnings(rm(".name", envir = env))
  return(res)
}


#' SGE qsub wrapper function
#'
#' @param X
#' @param FUN
#' @param keep_files
#' @param max_cores
#' @param jobname
#' @param tmpdir
#' @param env
#' @param ...
#'
#' @export
#'
#' @examples
qapply <- function(X, FUN, keep_files = FALSE, max_cores = 30, jobname = NULL, tmpdir = NULL, env = NULL, ...){

  ### Input ----

  Rscript <- system(paste0("source ~/.bashrc; which Rscript"), intern = TRUE)
  if (length(Rscript) == 0) stop("Error: 'Rscript' not found!")
  if (length(system(paste0("source ~/.bashrc; which qsub"), intern = TRUE)) == 0) stop("Error: 'qsub' not found!")
  if (is.null(jobname)) jobname <- "RPAR"

  jobnames <- system(paste0("source ~/.bashrc; qstat"), intern = TRUE)
  i <- 1
  jobname0 <- jobname
  while (jobname %in% jobnames){
    jobname <- paste0(jobname0, i)
    i <- i + 1
  }

  if (is.null(tmpdir)){
    tmpdir <- paste0(jobname, "_tmp")
    i <- 1
    while (dir.exists(tmpdir)){
      tmpdir <- paste0(paste0(jobname, "_tmp"), i)
      i <- i + 1
    }
  }

  if (!dir.exists(tmpdir) & keep_files == FALSE) on.exit(unlink(tmpdir, recursive = TRUE))
  if (!dir.exists(tmpdir)) dir.create(tmpdir)

  input <- list()
  input$FUN <- FUN
  input$list <- X
  input$args <- list(...)

  ncores <- min(length(input$list), max_cores)

  if (is.null(env)) env <- parent.frame()
  input$env <- env

  input$envs <- lapply(search(), as.environment)
  input$envs[[1]] <- as.environment(as.list(input$envs[[1]], all.names=TRUE)) # redefine global env
  input$envs <- rev(input$envs)

  saveRDS(object = input, file = file.path(tmpdir, "rscript_input.rds"))



  ### Write Rscript ----

  rscript <- file.path(tmpdir, "rfun.R")

  # HEADER
  cat(paste0("#!", Rscript), "\n", file = rscript)
  cat(paste0("#$ -S ", Rscript), "\n", file = rscript, append = TRUE)
  cat(paste0("#$ -pe smp ", ncores), "\n", file = rscript, append = TRUE)
  cat("#$ -cwd", "\n", file = rscript, append = TRUE)
  cat("#$ -V", "\n", file = rscript, append = TRUE)
  cat(paste0("#$ -N ", jobname), "\n", file = rscript, append = TRUE)
  cat(paste0("#$ -o ", file.path(tmpdir, "qsub.log")), "\n", file = rscript, append = TRUE)
  cat(paste0("#$ -e ", file.path(tmpdir, "qsub.err")), "\n", file = rscript, append = TRUE)


  # MAIN
  cat(paste0("input <- readRDS(file.path('", tmpdir, "', 'rscript_input.rds'))"), "\n", file = rscript, append = TRUE)
  cat(paste0('cl <- parallel::makeForkCluster(', ncores, ')'), "\n", file = rscript, append = TRUE)
  cat("on.exit(parallel::stopCluster(cl))", "\n", file = rscript, append = TRUE)
  cat("for (i in 1:length(input$envs)){ parallel::clusterExport(cl = cl, varlist = ls(input$envs[[i]]), envir = input$envs[[i]]) }",
      "\n", file = rscript, append = TRUE)
  cat("par_args <- c('cl' = list(cl), 'X' = list(input$list), 'fun' = input$FUN, input$args)", "\n", file = rscript, append = TRUE)
  cat("results <- do.call(parallel::parLapply, par_args)", "\n", file = rscript, append = TRUE)
  cat(paste0("saveRDS(object = results, file = file.path('", tmpdir, "', 'results.rds'))"), "\n", file = rscript, append = TRUE)


  ### Run Rscript ----
  system(paste0("source ~/.bashrc; qsub ", rscript), wait = TRUE)
  system(paste0("source ~/.bashrc; while qstat -s pr | grep -q -w ", jobname, "; do sleep 3; done"), wait = TRUE)


  ### Get results ----
  if (!file.exists(file.path(tmpdir, "results.rds"))){
    system(paste0("source ~/.bashrc; cat ", file.path(tmpdir, "qsub.err")))
    stop(paste0("Error: '", file.path(getwd(), tmpdir, 'results.rds'), "' not found."))
  } else {
    results <- readRDS(file.path(tmpdir, "results.rds"))
  }

  return(results)
}

