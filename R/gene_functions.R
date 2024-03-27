



#' Gene ID conversion
#'
#' @param ids
#' @param from
#' @param to
#' @param annotation
#' @param multiVals first, collapse
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' convertGeneIDs(c("3098", "3099"))
#' convertGeneIDs(c("3098", "3099"), to = "GENENAME")
#' convertGeneIDs(c("3098", "3099"), to = "GENETYPE")
convertGeneIDs <- function(ids, from = "ENTREZID", to = "SYMBOL", annotation = org.Hs.eg.db::org.Hs.eg.db, multiVals = "collapse", ...){

  stopifnot(requireNamespace("AnnotationDbi"))

  ktys <- AnnotationDbi::keytypes(annotation)
  if (!from %in% ktys | !to %in% ktys){
    print(paste0("Available ID types are: ", paste(ktys, collapse = ", ")))
    stop("Error: Keytype not found!")
  }

  if (multiVals == "collapse"){ multiVals <- function(x){paste0(x, collapse = "|")}}

  res <- AnnotationDbi::mapIds(x = annotation, keys = ids, column = to, keytype = from, multiVals = multiVals, ...)
  res[res == "NA"] <- NA
  res
}




#' Convert gene IDs between species
#'
#' @param genes
#' @param taxon
#' @param direction
#' @param id_type
#' @param na.omit
#'
#' @return
#' @export
#'
#' @examples
#' convertGeneSpecies("Egfr", taxon = "10090")
#' convertGeneSpecies("EGFR", taxon = "10090", direction = "to")
convertGeneSpecies <- function(genes, taxon = "10090", direction = "from", id_type = "symbol", na.omit = TRUE){

  stopifnot(requireNamespace("babelgene"))

  id_type <- tolower(id_type)

  tax_avail <- babelgene:::orthologs_df$taxon_id |> unique()
  if (!taxon %in% tax_avail){
    message("Available taxon IDs:")
    message(paste(tax_avail, collapse = ", "))
    stop("Error: Taxon ID not found.")
  }

  id_avail <- c("symbol", "entrez", "ensembl")
  if (!id_type %in% id_avail){
    message("Available gene ID types:")
    message(paste(id_avail, collapse = ", "))
    stop("Error: Gene ID type not found.")
  }

  df <- babelgene:::orthologs_df |> subset(taxon_id == taxon)

  if (direction == "from"){

    genes <- df[match(genes, df[[id_type]]),][[paste0("human_", id_type)]]

  } else if (direction == "to"){
    genes <- df[match(genes, df[[paste0("human_", id_type)]]),][[id_type]]
  }

  if (na.omit == TRUE){
    genes <- genes |> na.omit()
    genes <- genes |> as.character()
  }

  genes
}




#' Get gene sets from msigdbr
#'
#' @param collections Collections (msigdbr::msigdbr_collections()) and subcollections
#' @param species
#' @param id_type
#' @param format
#'
#' @return
#' @export
#'
#' @examples getGeneSets(c("H", "C2|CP:KEGG|CP:REACTOME"))
getMSigDB <- function(collections = c("H", "C2|CP:KEGG"), species = "human", id_type = "gene_symbol", format = c("list", "dataframe", "GeneSetCollection"), ...){

  stopifnot(requireNamespace("msigdbr"))
  # msigdbr::msigdbr_species()
  # msigdbr::msigdbr_collections()

  if (!is.null(collections)){

    collections <- strsplit(collections, split = "|", fixed = TRUE)
    categories <- sapply(collections, function(x) x[1] )
    subcategories <- lapply(setNames(collections, categories), function(x) unique(x[-1]) )
    subcategories[sapply(subcategories, length) == 0] <- NA
    df <- stack(subcategories)[,c(2,1)]
    colnames(df) <- c("category", "subcategory")


    gslist <- lapply(1:nrow(df), function(i){
      categ <- as.character(df[i,1])
      subcateg <- df[i,2]
      if (is.na(subcateg)) subcateg <- NULL
      msigdbr::msigdbr(species = species, category = categ, subcategory = subcateg)
    })

    genesets <- Reduce(x = gslist, f = rbind)

  } else {

    genesets <- msigdbr::msigdbr(species = species)
  }


  df <- unique(data.frame(genesets[,c(id_type, "gs_name")]))
  colnames(df) <- c("gene", "term")
  GS <- convertGeneSets(df, from = "dataframe", to = format, ...)

  GS
}





#' Get gene/metabolite sets from metabolicatlas
#'
#' @param type
#'
#' @return
#' @export
#'
#' @examples
getMetabolicAtlas <- function(type = "genes"){

  stopifnot(tolower(type) %in% c("genes", "metabolites"))

  json_subsystems <- rjson::fromJSON(paste(readLines("https://metabolicatlas.org/api/v2/maps/listing?model=HumanGem"), collapse=""))
  subsystems <- sapply(json_subsystems$subsystems, function(tmp) tmp$id )
  subsystem_data <- lapply(setNames(subsystems, subsystems), function(tmp){
    url <- paste0("https://metabolicatlas.org/api/v2/subsystems/", tmp,"?model=HumanGem")
    rjson::fromJSON(paste(readLines(url, warn = FALSE), collapse=""))
  })


  if (type == "genes"){
    HUMAN1 <- sapply(subsystem_data, function(tmp){
      unique(unlist( sapply(tmp$genes, function(tmp2){ tmp2$name }) ))
    })
  }

  if (type == "metabolites"){
    HUMAN1 <- sapply(subsystem_data, function(tmp){
      unique(unlist( sapply(tmp$metabolites, function(tmp2){ tmp2$id }) ))
    })
  }

  names(HUMAN1) <- paste0("HUMAN1_", names(HUMAN1))
  HUMAN1
}




#' Get gene sets from MitoCarta3
#'
#' @param file
#'
#' @return
#' @export
#'
#' @examples
getMitoCarta <- function(file = "Human.MitoCarta3.0.xls", ...){

  db <- "ftp://ftp.broadinstitute.org/distribution/metabolic/papers/Pagliarini/MitoCarta3.0/Human.MitoCarta3.0.xls"
  download.file(db, destfile = file)
  mitocarta <- as.data.frame(readxl::read_xls(file, sheet = 2))
  gs <- mitocarta$MitoCarta3.0_MitoPathways %>% strsplit(split = ">|\\|") %>% sapply(trimws) %>% sapply(tolower)
  gsnames <- gs %>% unlist() %>% as.character() %>% unique()
  gsnames <- gsnames[!is.na(gsnames)]
  gsnames <- gsnames[gsnames != "0"]
  mcgs <- lapply(gsnames,function(x) mitocarta$Symbol[ sapply(gs,function(y) x %in% y )] )
  mcgs <- mcgs[!is.na(gsnames)]
  names(mcgs) <- gsnames[!is.na(gsnames) & nchar(gsnames) > 1]
  mcgs <- stack(mcgs)[,c(2,1)]
  colnames(mcgs) <- c("term", "gene")
  mcgs

}








#' Get GO gene sets
#'
#' @param db
#' @param keytype
#' @param ontology
#' @param minsize
#'
#' @return
#' @export
#'
#' @examples
getGOgenes <- function(db = org.Hs.eg.db::org.Hs.eg.db, keytype = "SYMBOL", ontology = c("BP", "MF", "CC"), minsize = 3){

  go2ont <- AnnotationDbi::Ontology(GO.db::GOTERM)

  godf <- stack(AnnotationDbi::mapIds(db, keys = unique(names(go2ont)), column = keytype, keytype = "GOALL", multiVals = "list"))
  colnames(godf) <- c("gene", "id")
  godf$gene <- as.character(godf$gene)
  godf$id <- as.character(godf$id)
  godf <- subset(godf, !is.na(gene) & !is.na(id))

  go2term <- AnnotationDbi::Term(GO.db::GOTERM[unique(godf$id)])

  godf$term <- go2term[godf$id]
  godf$ont <- go2ont[godf$id]

  godf <- subset(godf, ont != "universal")
  godf <- subset(godf, ont %in% ontology)

  godf$size <- table(godf$id)[godf$id]
  godf <- subset(godf, size >= minsize)

  data.frame(godf[,c("gene", "term", "id", "ont", "size")], row.names = NULL)
}





#' Convert gene sets between different formats
#'
#' @param genesets
#' @param from
#' @param to
#' @param term
#' @param gene
#'
#' @return
#' @export
#'
#' @examples
convertGeneSets <- function(genesets, from = NULL, to = "list", term = term, gene = gene){

  stopifnot(requireNamespace("GSEABase"))

  term <- rlang::enquo(term)
  gene <- rlang::enquo(gene)

  if (is.null(from)){
    if ("data.frame" %in% class(genesets)) from <- "dataframe"
    if ("list" %in% class(genesets)) from <- "list"
    if ("GeneSetCollection" %in% class(genesets)) from <- "GeneSetCollection"
  }

  from <- trimws(tolower(from[1]))
  to <- trimws(tolower(to[1]))

  gsc <- FALSE
  if (to %in% c("list", "l")) to <- "list"
  if (to %in% c("dataframe", "df", "data.frame", "table")) to <- "dataframe"
  if (to %in% "genesetcollection"){
    to <- "list"
    gsc <- TRUE
  }

  if (from == "GeneSetCollection" & to == "GeneSetCollection") return(genesets)


  # from list to dataframe
  if (from == "list" & to == "dataframe"){
    genesets <- utils::stack(genesets)[,c(2,1)]
    colnames(genesets) <- c(rlang::as_name(term), rlang::as_name(gene))
  }

  # from dataframe to list
  if (from == "dataframe" & to == "list"){
    genesets <- dplyr::select(genesets, !!gene, !!term)
    genesets <- unstack(unique(genesets))
  }

  # from list to GeneSetCollection
  if (gsc == TRUE){
    genesets <- lapply(setNames(seq_along(genesets), names(genesets)), function(i){
      GSEABase::GeneSet(genesets[[i]], setName = names(genesets[i]))
    })
    genesets <- GSEABase::GeneSetCollection(genesets)
  }


  genesets
}




#' Get random genes
#'
#' @param n
#' @param keytype
#' @param dbi
#' @param replace
#' @param prob
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' rgenes(100)
rgenes <- function(n = 20, keytype = "SYMBOL", dbi = org.Hs.eg.db::org.Hs.eg.db, replace = FALSE, prob = NULL, ...){
  AnnotationDbi::keys(dbi, keytype = keytype, ...) |> sample(size = n, replace = replace, prob = prob)
}




#' Read GMT files
#'
#' @param file
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
readGMT <- function(file, ...){

  stopifnot(requireNamespace("GSEABase"))

  gmt <- GSEABase::getGmt(file, ...)
  GSEABase::geneIds(gmt)

}



#' Write GMT files
#'
#' @param genesets
#' @param file
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
writeGMT <- function(genesets, file, ...){

  stopifnot(requireNamespace("GSEABase"))

  if ("genesetcollection" %in% tolower(class(genesets))){
    genesets <- convertGeneSets(genesets, to = "GeneSetCollection")
  }

  GSEABase::toGmt(genesets, file, ...)

}


#' Read gtf files
#'
#' @param file
#' @param columns
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
readGTF <- function(file, columns = NULL, ...){

  gtf <- rtracklayer::import(file)
  gtf <- subset(gtf, ...)

  if (is.null(columns)){
    columns <- c("gene_name", "seqnames", "gene_type")
    cat(paste0("Returning columns: ", paste(columns, collapse = ", ")))
    cat(paste0("Available columns: ", paste(colnames(as.data.frame(head(gtf))), collapse = ", ")))
  } else if ("all" %in% tolower(columns)){
    columns <- colnames(as.data.frame(head(gtf)))
  }

  df <- as.data.frame(gtf)[,columns, drop = FALSE]
  df <- dplyr::distinct(df)
  df
}













#' Retrieve data from BiomaRt
#'
#' @param genes
#' @param idtype
#' @param attributes
#' @param biotype
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
getBiomaRt <- function(genes, idtype = "hgnc_symbol", attributes = "transcript_biotype", biotype = "protein_coding", ...){

  stopifnot(requireNamespace("biomaRt", quietly = TRUE))

  mart <- biomaRt::useMart(biomart = "ensembl", "hsapiens_gene_ensembl")
  results <- biomaRt::getBM(attributes = c(idtype, attributes), filters = idtype, values = genes, mart = mart, ...)

  results
}









