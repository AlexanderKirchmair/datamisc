




#' Load Seurat datasets
#'
#' @param ds Dataset name
#'
#' @return
#' @export
#'
#' @examples
#' get_seu_dataset()
get_seu_dataset <- function(ds = "pbmc3k"){
  stopifnot(requireNamespace("SeuratData"))
  if (!ds %in% SeuratData::InstalledData()[1]){
    SeuratData::InstallData(ds)
  }
  SeuratData::LoadData(ds)
}



#' Downsample Seurat object
#'
#' @param seu Seurat object
#' @param group_by Cluster
#' @param n_cells Number of target cell
#' @param min_cells Minimum number of cells per group
#' @param seed Random seed
#'
#' @return
#' @export
#'
#' @examples
#' get_seu_dataset() |> seu_downsample(group_by = "seurat_annotations", n_cells = 100, min_cells = 10)
seu_downsample <- function(seu, group_by = "celltype", n_cells = NULL, min_cells = NULL, seed = 123){

  if (!is.null(seed)) set.seed(seed)

  if (is.null(n_cells)) n_cells <- round(ncol(seu) * 0.1, digits = -2)
  if (is.null(min_cells)) min_cells <- round(n_cells / length(unique(seu@meta.data[[group_by]])) * 0.3)

  if (sum(is.na(seu@meta.data[[group_by]])) > 0){
    warning("Missing values present in 'group_by', converting to factor")
    seu@meta.data[[group_by]] <- add_level(seu@meta.data[[group_by]], is.na, "NA")
  }

  groups_orig <- table(seu@meta.data[[group_by]]) |> sort(decreasing = TRUE)
  groups <- groups_orig

  if (length(groups) * min_cells > n_cells){
    stop("The total number of cells 'n_cells' must be larger than 'min_cells' * number of groups!")
  }

  print(paste0("Downsampling cells with n_cells=", n_cells, " and min_cells=", min_cells, ""))

  n_cells_per_group <- groups
  while (!all(n_cells_per_group >= min_cells) | sum(groups) > n_cells){
    groups_to_downsample <- groups[groups > min_cells]
    groups_not_downsample <- groups[groups <= min_cells]
    n_cells_tagret <- n_cells - sum(groups_not_downsample)
    frac_target <- n_cells_tagret / sum(groups_to_downsample)
    n_cells_per_group <- round(groups_to_downsample * frac_target)

    while (sum(n_cells_per_group) != n_cells_tagret){
      if (sum(n_cells_per_group) > n_cells_tagret){
        n_cells_per_group[1] <- n_cells_per_group[1] - 1
      } else {
        n_cells_per_group[1] <- n_cells_per_group[1] + 1
      }
    }

    n_cells_per_group[n_cells_per_group < min_cells] <- min_cells
    groups <- c(n_cells_per_group, groups_not_downsample)[names(groups)]
  }

  cells <- lapply(setNames(names(groups), names(groups)), function(g){
    sample(which(seu@meta.data[[group_by]] == g), groups[g])
  })

  seu_small <- subset(seu, cells = sort(unique(unlist(cells))))

  groups_small <- table(seu_small@meta.data[[group_by]]) |> sort(decreasing = TRUE)
  stopifnot(sum(groups_small) <= n_cells)
  stopifnot(all(names(groups_small) %in% names(groups_orig)))

  groups_small <- groups_small[!groups_small %in% groups_small[duplicated(groups_small)]]
  if (!all(names(groups_small) == names(groups_orig)[names(groups_orig) %in% names(groups_small)])){
    warning("Warning: Relative group sizes are not preserved!")
  }

  return(seu_small)
}






#' Find markers for groups in a Seurat object
#'
#' @param seu
#' @param group_by
#' @param group_oi
#' @param group_ref
#' @param assay
#' @param slot
#' @param logfc_thresh
#' @param min_pct
#' @param padj_thresh
#' @param test_use
#' @param rename
#' @param n_cores
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' get_seu_dataset() |> find_markers(group_by = "seurat_annotations", logfc_thresh = 1, min_pct = 0.5, padj_thresh = 0.05)
#' get_seu_dataset() |> find_markers(group_by = "seurat_annotations", group_oi = "DC", group_ref = "NK")
find_markers <- function(seu, group_by = NULL, group_oi = NULL, group_ref = NULL, assay = NULL, slot = "data", logfc_thresh = 0.1, min_pct = 0.1, padj_thresh = 0.05, test_use = "wilcox", rename = TRUE,  n_cores = 1, ...){

  stopifnot(requireNamespace("Seurat"))

  if (!is.null(group_by)){
    Idents(seu) <- seu@meta.data[[group_by]]
  }

  if (n_cores > 1){
    future::plan("multisession", workers = n_cores)
    on.exit(future::plan("sequential"))
  }

  if (is.null(group_oi) & is.null(group_ref)){
    de_markers <- FindAllMarkers(object = seu, assay = assay, logfc.threshold = logfc_thresh, test.use = test_use, slot = slot, min.pct = min_pct, min.diff.pct = -Inf, ...)

  } else {
    de_markers <- FindMarkers(object = seu, ident.1 = group_oi, ident.2 = group_ref, assay = assay, logfc.threshold = logfc_thresh, test.use = test_use, slot = slot, min.pct = min_pct, min.diff.pct = -Inf, ...)
    de_markers$gene <- rownames(de_markers)

    if (is.null(group_ref)){
      de_markers$group <- group_oi
    } else {
      de_markers$group <- paste0(group_oi, "_vs_", group_ref)
    }

  }

  if (all(c("pct.1", "pct.2") %in% colnames(de_markers))) {
    de_markers$pct.diff <- de_markers$pct.1 - de_markers$pct.2
  }

  de_markers$cluster <- as.character(de_markers$cluster)

  de_markers <- de_markers |> dplyr::relocate(cluster, gene)
  de_markers <- de_markers |> dplyr::relocate(p_val_adj, .after = p_val)

  if (all(c("p_val", "avg_log2FC", "pct.diff") %in% colnames(de_markers))){
    de_markers <- de_markers |> dplyr::arrange(p_val, desc(abs(avg_log2FC)), desc(pct.diff))
  }

  if (!is.null(padj_thresh)){
    de_markers <- de_markers |> dplyr::filter(p_val_adj <= padj_thresh)
  }

  if (rename == TRUE){
    de_markers <- de_markers |> dplyr::rename(group = cluster, pval = p_val, padj = p_val_adj, log2FC = avg_log2FC, pct1 = pct.1, pct2 = pct.2, pctdiff = pct.diff)
  }

  rownames(de_markers) <- NULL
  de_markers
}




#' Extract top markers from find_markers() results
#'
#' @param markers Dataframe
#' @param n Number of top features
#' @param per_group Get features per group
#' @param distinct Count each feature only once
#' @param group_col
#' @param feature_col
#' @param filter_pattern
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' get_seu_dataset() |> find_markers(group_by = "seurat_annotations") |> get_top_markers()
get_top_markers <- function(markers, n = 5, per_group = TRUE, distinct = FALSE, group_col = "group", feature_col = "gene", filter_pattern = NULL, ...){

  stopifnot(!is.null(markers[[feature_col]]))
  stopifnot(!is.null(markers[[group_col]]))

  markers[[group_col]] <- as.character(markers[[group_col]])

  # sort
  if (length(enquos(...)) != 0){
    markers <- markers |> dplyr::arrange(...)
  }

  # filter
  markers <- markers |> dplyr::select(!!group_col, !!feature_col)
  if (!is.null(filter_pattern)){
    markers <- markers[!grepl(pattern = filter_pattern, x = markers[[feature_col]]),]
  }
  if (distinct == TRUE){
    markers <- markers[!duplicated(markers[[feature_col]]),,drop=FALSE]
  }

  # split
  if (per_group == TRUE){
    markers <- markers |> split(f = markers[[group_col]])
    top_markers <- lapply(markers, function(df) df |> dplyr::slice_head(n = n) |> dplyr::pull(!!feature_col))
  } else {
    top_markers <- markers |> dplyr::slice_head(n = n) |> dplyr::pull(!!feature_col)
  }

  top_markers
}



#' DotPlot wrapper for Seurat
#'
#' @param seu
#' @param group_by
#' @param features
#' @param dot_scale
#' @param col_mid
#' @param scale
#' @param col_low
#' @param col_high
#' @param col_min
#' @param col_max
#' @param fontsize
#' @param lwd
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' get_seu_dataset() |> seu_dotplot(group_by = "seurat_annotations", features = c( "S100A9", "S100A8", "LGALS2", "CD79A", "GZMB", "HLA-DQA1", "CCL5", "NKG7"))
seu_dotplot <- function(seu, group_by = NULL, features = NULL, dot_scale = NULL, col_mid = NULL, scale = FALSE, col_low = "grey98", col_high = "#0000b5", col_min = NA, col_max = NA, fontsize = 12, lwd = 0.5, ...){

  stopifnot(requireNamespace("Seurat"))

  if (sum(is.na(seu@meta.data[[group_by]])) > 0){
    warning("Missing values present in 'group_by', converting to factor")
    seu@meta.data[[group_by]] <- add_level(seu@meta.data[[group_by]], is.na, "NA")
  }

  if (is.null(dot_scale)){
    dot_scale <- 90 / sqrt(length(unique(seu@meta.data[[group_by]])) * length(features))
  }

  gg <- Seurat::DotPlot(seu, group.by = group_by, cols = c("blue", "grey97", "red"), features = features, dot.scale = dot_scale, scale = scale, ...) +
    ggplot2::theme(text = element_text(size = fontsize),
                   axis.text.x = element_text(size = fontsize, angle = 90, hjust = 1, vjust = 0.5),
                   axis.text.y = element_text(size = fontsize),
                   legend.box = "horizontal")

  gg <- gg + xlab("") + ylab("")

  if (is.null(col_mid)){
    gg <- gg + scale_color_gradient(low = col_low, high = col_high, limits = c(col_min, col_max), oob = scales::squish)
  } else {
    gg <- gg + scale_color_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = 0, oob = scales::squish)
  }

  gg <- gg + theme(axis.line = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7))

  gg <- gg + guides(color = guide_colorbar(title = "mean", frame.colour = "black", ticks.colour = "black", frame.linewidth = lwd),
                    size = guide_legend(title = "cells (%)"))
  gg
}






#' UMAP plot wrapper for Seurat
#'
#' @param seu
#' @param features
#' @param groups
#' @param label
#' @param alpha
#' @param assay
#' @param slot
#' @param colors
#' @param min_cutoff
#' @param max_cutoff
#' @param reduction
#' @param pt_shape
#' @param pt_size
#' @param title_size
#' @param axis_size
#' @param axis_text_size
#' @param legend_key_size
#' @param legend_text_size
#' @param lwd
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' seu <- get_seu_dataset() |> ScaleData() |> FindVariableFeatures() |> RunPCA() |> RunUMAP(dims = 1:10)
#' seu_umap(seu, features = "S100A9")
#' seu_umap(seu, groups = "seurat_annotations")
seu_umap <- function(seu, features = NULL, groups = NULL, label = FALSE, alpha = 1, assay = NULL, slot = NULL, colors = NULL, min_cutoff = "q05", max_cutoff = "q95", reduction = "umap",
                     pt_shape = 16, pt_size = 0.5, title_size = 18, axis_size = 16, axis_text_size = 14, legend_key_size = 18, legend_text_size = 18, lwd = 0.5, ...){

  stopifnot(requireNamespace("Seurat"))

  update_gg <- function(ggobj){
    ggobj$layers[[1]]$aes_params[["shape"]] <- pt_shape
    ggobj <- ggobj + xlab("UMAP1") + ylab("UMAP2")
    if (!is.null(groups)){
      ggobj <- ggobj + guides(color = guide_legend(override.aes = list(size = legend_key_size/ggplot2:::.pt)))
    }
    ggobj
  }

  if (!is.null(assay)){
    DefaultAssay(seu) <- assay
  }

  if (!is.null(groups)){
    gg <- Seurat::DimPlot(seu, group.by = groups, cols = colors,
                  pt.size = pt_size, alpha = alpha,
                  reduction = reduction, label = label, ...)
  }

  if (!is.null(features)){
    if (is.null(colors)){
      colors <- c("#c9c9c9", "#0000ff")
    }
    gg <- Seurat::FeaturePlot(seu, features = features, cols = colors,
                      pt.size = pt_size, alpha = alpha, order = TRUE,
                      slot = slot, min.cutoff = min_cutoff, max.cutoff = max_cutoff, reduction = reduction, ...)
    gg <- gg + guides(color = guide_colorbar(title = "", frame.colour = "black", ticks.colour = "black", frame.linewidth = lwd))
  }

  for (i in 1:length(gg)){
    gg[[i]] <- update_gg(gg[[i]])
  }

  gg & theme(axis.title = element_text(size = axis_size),
             axis.text = element_text(size = axis_text_size),
             plot.title = element_text(size = title_size),
             legend.text = element_text(size = legend_text_size))
}



