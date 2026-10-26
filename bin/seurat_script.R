#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(Seurat)
    library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 1 && args[[1]] %in% c("--help", "-h")) {
    cat("Usage: Rscript seurat_script.R --input <filtered_feature_bc_matrix> --output <output_prefix>\n")
    quit(status = 0)
}

get_arg <- function(flag) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) {
        stop(sprintf("Missing required argument: %s", flag))
    }
    args[[i + 1]]
}

input_dir <- get_arg("--input")
output_prefix <- get_arg("--output")

if (!dir.exists(input_dir)) {
    stop(sprintf("Input matrix directory does not exist: %s", input_dir))
}

options(future.globals.maxSize = 8 * 1024^3)

counts <- Read10X(data.dir = input_dir)
if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
    } else {
        counts <- counts[[1]]
    }
}

obj <- CreateSeuratObject(
    counts = counts,
    project = basename(output_prefix),
    min.cells = 3,
    min.features = 200
)

obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

# Conservative generic QC defaults. Adjust for the project if needed.
obj <- subset(
    obj,
    subset = nFeature_RNA >= 200 &
             nFeature_RNA <= 7500 &
             percent.mt <= 20
)

if (ncol(obj) < 10) {
    stop(sprintf("Too few cells remain after QC for %s: %d", output_prefix, ncol(obj)))
}

obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(2000, nrow(obj)))
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)

n_pcs <- min(30, ncol(obj) - 1, length(VariableFeatures(obj)) - 1)
if (n_pcs < 2) {
    stop(sprintf("Too few dimensions remain for PCA in %s", output_prefix))
}

obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = n_pcs, verbose = FALSE)
dims_use <- seq_len(n_pcs)
obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, dims = dims_use, verbose = FALSE)

pdf(paste0(output_prefix, ".pdf"), width = 11, height = 8.5)
print(VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
print(ElbowPlot(obj, ndims = n_pcs))
print(DimPlot(obj, reduction = "pca", group.by = "seurat_clusters"))
print(DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE))
print(DimPlot(obj, reduction = "umap", group.by = "percent.mt"))
dev.off()

saveRDS(obj, file = paste0(output_prefix, ".rds"))

cat(sprintf("Seurat analysis complete: %s\n", output_prefix))
cat(sprintf("Cells retained: %d\n", ncol(obj)))
cat(sprintf("Genes retained: %d\n", nrow(obj)))
