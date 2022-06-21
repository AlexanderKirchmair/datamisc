


df2graph <- function(edges, nodes, directed = TRUE){

  # add checks

  graph <- igraph::graph_from_data_frame(d = edges, vertices = nodes, directed = directed)
  graph <- tidygraph::as_tbl_graph(graph, stringsAsFactors = FALSE)
  graph
}


node2edge <- function(graph, col = "stat", direction = "out", FUN = NULL, ...){

  edf <- edgeData(graph)
  ndf <- nodeData(graph)

  val <- setNames(ndf[[col]], ndf$name)

  if (direction == "out"){
    edf[[col]] <- val[edf$from]
  } else if (direction == "in"){
    edf[[col]] <- val[edf$to]
  } else if (direction == "both"){
    if (is.null(FUN)) stop("Provide FUN to combine values!")
    edf[[col]] <- apply(cbind(val[edf$from], val[edf$to]), 1, FUN, ...)
  } else {
    warning("'Direction' muste be one of: out, in. both")
  }

  df2graph(edf, ndf)
}




# map node data to edge source nodes
sources <- function(graph, ...){

  edf <- edgeData(graph, ...)
  ndf <- nodeData(graph)

  if (is.null(edf$from) | is.null(ndf$name)){ stop("Error: missing node/edge ids!") }
  df <- ndf[match(edf$from, ndf$name),, drop = FALSE]

  rownames(df) <- rownames(edf)
  stopifnot(nrow(edf) == nrow(df))

  df
}



targets <- function(graph, ...){

  edf <- edgeData(graph, ...)
  ndf <- nodeData(graph)

  if (is.null(edf$to) | is.null(ndf$name)){ stop("Error: missing node/edge ids!") }
  df <- ndf[match(edf$to, ndf$name),, drop = FALSE]

  rownames(df) <- rownames(edf)
  stopifnot(nrow(edf) == nrow(df))

  df
}









edgeData <- function(graph, ...){
  edf <- igraph::as_data_frame(x = graph, what = "edges")
  edf <- subset(edf, ...)
  edf
}


nodeData <- function(graph, ...){
  ndf <- igraph::as_data_frame(x = graph, what = "vertices")
  ndf <- subset(ndf, ...)
  ndf
}



lcc <- function(graph){

  # Function that returns the largest connected component of a network

  subgraphs <- igraph::decompose.graph(graph) # decompose graph into disconnected subgraphs/components
  subgraphs_nodes <- sapply(subgraphs, igraph::vcount) # get the number of nodes for each subgraph

  if (any(subgraphs_nodes)){

    ixlargest <- (1:length(subgraphs_nodes))[subgraphs_nodes == max(subgraphs_nodes)]
    res <- subgraphs[[ixlargest[1]]] # select largest component

    # messages
    Vkept <- length(igraph::as_ids(igraph::V(res)))
    Vorig <- length(igraph::as_ids(igraph::V(res)))


    message(paste("Number of disconnected components is ", length(subgraphs), ".", sep = ""))
    message(paste("The largest component has ", Vkept, " nodes.", sep = ""))
    message(paste(Vorig - Vkept, " nodes were removed.", sep = ""))


  } else {
    print("No nodes found in graph.")
    res <- graph
  }


  res
}



rgraph <- function(n = 20){

  rnodes <- paste0(LETTERS, rep(1:3, each = length(LETTERS)))[1:n]
  graph <- tidygraph::play_erdos_renyi(n = n, p = 0.1, loops = TRUE)
  edf <- igraph::as_data_frame(graph)
  edf <- rbind(edf, edf[sample(1:nrow(edf), round(nrow(edf)*0.3)),])
  graph <- tidygraph::as_tbl_graph(igraph::graph_from_data_frame(edf))

  graph
}










diffuse <- function(graph, nstat, estat, na = 0, maxiter = 1000, tresh = 10^-6, decay = 1, fixed = TRUE, FUN = NULL, ...){

  ndf <- nodeData(graph)
  edf <- edgeData(graph)

  ERR <- function(v, vdiff){
    v <- as.numeric(v)
    mean(abs( (v - as.numeric(vdiff))/v), na.rm = TRUE, trim = 0.1)
  }



  # Node diffusion ---

  v <- ndf[[nstat]]
  v_fixed <- !is.na(v)
  v_free <- is.na(v)
  v[is.na(v)] <- na

  A <- igraph::as_adjacency_matrix(graph, sparse = TRUE) * decay
  n <- sparseMatrixStats::rowSums2(A != 0, na.rm = TRUE)

  if (is.null(FUN)){
    cat(crayon::blue("Network diffusion by matrix multiplication...\n"))

    i <- 0
    err <- Inf

    if (fixed == TRUE){
      while (i < maxiter & err > tresh){
        i <- i + 1
        vdiff <- (A %*% v)/n
        err <- ERR(v[v_free], vdiff[v_free])
        v[v_free] <- vdiff[v_free]
      }
    } else {
      while (i < maxiter & err > tresh){
        i <- i + 1
        vdiff <- (A %*% v)/n
        err <- ERR(v, vdiff)
        v <- vdiff
      }
    }


    cat(crayon::blue(paste0("Stopped after ", i, " iterations, error = ", signif(err, 5))))

  } else {
    cat(crayon::blue("Network diffusion using custom function...\n"))


      vdiff <- FUN(v, A, ...)

    cat(crayon::blue(paste0("Stopped after ", i, " iterations, error = ", signif(err, 5))))

  }

  ndf[[nstat]] <- as.numeric(v)
  df2graph(edf, ndf)
}







reverse_directions <- function(graph){

  from <- edgeData(graph)$from
  to <- edgeData(graph)$to

  edgeData(graph)$from <- to
  edgeData(graph)$to <- from

  graph
}







subset.igraph <- function(graph, type, ...){

  ndf <- nodeData(graph)
  edf <- edgeData(graph)

  if (type == "nodes"){
    ndf <- subset(ndf, ...)
    edf <- subset(edf, from %in% ndf$name & to  %in% ndf$name)

  } else if (type == "edges"){
    edf <- subset(edf, ...)
    ndf <- subset(ndf, name %in% edf$from | name %in% edf$to)
  }

  df2graph(edf, ndf, directed = igraph::is.directed(graph))
}





layout_with_ggrepel <- function(graph, labels = NULL, init = NULL, scale = TRUE, point_size = 0, point_padding_x = 0.02, point_padding_y = 0.01, max_iter = 10^5, max_time = 30, ...){

  stopifnot(requireNamespace("ggrepel"))

  if (is.null(labels)) labels <- ""
  if (is.null(init)) init <- graphlayouts::layout_with_stress(graph)
  if (length(point_size) == 1) point_size <- rep(point_size, nrow(init))

  xy <- init
  if (scale == TRUE) xy <- apply(xy, 2, rangescale)
  colnames(xy) <- c("x", "y")

  xybox <- cbind(xy[,"x"] - point_padding_x/2 , xy[,"y"] - point_padding_x/2,
                 xy[,"x"] + point_padding_x/2, xy[,"y"] + point_padding_y/2)

  xyrep <- ggrepel:::repel_boxes2(data_points = xy,
                                  point_size = point_size,
                                  point_padding_x = point_padding_x,
                                  point_padding_y = point_padding_y,
                                  boxes = xybox,
                                  xlim = range(c(xybox[,1], xybox[3])),
                                  ylim = range(c(xybox[,2], xybox[4])),
                                  hjust = rep(0.5, nrow(xy)),
                                  force_push = 10^-6,
                                  force_pull = 1.1*10^-6,
                                  vjust = rep(0.5, nrow(xy)),
                                  max_overlaps = Inf,
                                  max_time = max_time,
                                  max_iter = max_iter, ...)


  res <- xyrep[,c("x", "y")]
  res <- as.matrix(res)
  if (scale == TRUE) res <- apply(res, 2, rangescale)
  dimnames(res) <- dimnames(xy)
  res

}






leafs <- function(graph, mode = "all", multiple = FALSE, loops = FALSE, mindegree = 1, names = FALSE, logical = FALSE, ...){

  if (!multiple) graph <- igraph::simplify(graph, remove.multiple = !multiple, remove.loops = !loops, edge.attr.comb = "ignore")
  is_leaf <- igraph::degree(graph, mode = mode, loops = loops, ...) == mindegree

  if (names == TRUE){
    val <- V(graph)[ix]$name
  } else if (logical == TRUE){
    val <- is_leaf
  } else {
    val <- which(is_leaf)
  }

  unname(val)
}










