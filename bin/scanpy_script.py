#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(description="Standard Scanpy analysis of a Cell Ranger matrix")
    parser.add_argument("--input", required=True, help="Cell Ranger filtered_feature_bc_matrix directory")
    parser.add_argument("--output", required=True, help="Output prefix")
    return parser.parse_args()


def main():
    args = parse_args()
    input_dir = Path(args.input)
    output_prefix = Path(args.output)

    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input matrix directory does not exist: {input_dir}")

    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    adata = sc.read_10x_mtx(
        input_dir,
        var_names="gene_symbols",
        make_unique=True,
        cache=False,
    )

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Conservative generic QC defaults. Adjust for the project if needed.
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["pct_counts_mt"] <= 20, :].copy()
    adata = adata[adata.obs["n_genes_by_counts"] <= 7500, :].copy()

    if adata.n_obs < 10:
        raise RuntimeError(f"Too few cells remain after QC: {adata.n_obs}")

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata

    n_hvg = min(2000, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat")
    if adata.var["highly_variable"].sum() < 2:
        raise RuntimeError("Fewer than two highly variable genes were identified")

    sc.pp.scale(adata, max_value=10)
    n_pcs = min(30, adata.n_obs - 1, int(adata.var["highly_variable"].sum()) - 1)
    if n_pcs < 2:
        raise RuntimeError(f"Too few dimensions remain for PCA: {n_pcs}")

    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1), n_pcs=n_pcs)
    sc.tl.umap(adata)

    # Leiden is optional because leidenalg may not be installed in every Scanpy environment.
    try:
        sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2)
        cluster_key = "leiden"
    except Exception:
        cluster_key = None

    output_h5ad = output_prefix.with_suffix(".h5ad")
    output_pdf = output_prefix.with_suffix(".pdf")
    adata.write_h5ad(output_h5ad)

    with PdfPages(output_pdf) as pdf:
        fig, axes = plt.subplots(1, 3, figsize=(14, 4))
        axes[0].hist(adata.obs["n_genes_by_counts"], bins=50)
        axes[0].set_title("Genes per cell")
        axes[1].hist(adata.obs["total_counts"], bins=50)
        axes[1].set_title("Counts per cell")
        axes[2].hist(adata.obs["pct_counts_mt"], bins=50)
        axes[2].set_title("Mitochondrial percentage")
        for ax in axes:
            ax.set_xlabel("Value")
            ax.set_ylabel("Cells")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        axes[0].scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
                        c=adata.obs["n_genes_by_counts"], s=3, cmap="viridis")
        axes[0].set_title("UMAP: genes per cell")
        axes[1].scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
                        c=adata.obs["pct_counts_mt"], s=3, cmap="magma")
        axes[1].set_title("UMAP: mitochondrial percentage")
        for ax in axes:
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        if cluster_key is not None:
            codes = adata.obs[cluster_key].astype("category").cat.codes.to_numpy()
            fig, ax = plt.subplots(figsize=(7, 6))
            scatter = ax.scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
                                 c=codes, s=3, cmap="tab20")
            ax.set_title(f"UMAP: {cluster_key} clusters")
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
            pdf.savefig(fig)
            plt.close(fig)

    print(f"Scanpy analysis complete: {output_prefix}")
    print(f"Cells retained: {adata.n_obs}")
    print(f"Genes retained: {adata.n_vars}")
    print(f"H5AD output: {output_h5ad}")


if __name__ == "__main__":
    main()
