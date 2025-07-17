import scanpy as sc
import numpy as np

def process_adata(adata: sc.AnnData) -> None:
    """
    Process the an AnnData object by filtering cells and genes, computing graphs, and clustering.

    Parameters:
    adata (sc.AnnData): The AnnData object to process.
    """
    adata.var_names_make_unique()  # Ensure unique gene names
    adata.obs_names_make_unique()  # Ensure unique cell names
    
    # Filter cells and genes
    sc.pp.filter_cells(
        adata, min_genes=20
    )  # From (https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html)
    sc.pp.filter_genes(adata, min_cells=np.exp(4))  # From Roshan's workshop

    # Compute graphs
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=4000,
        batch_key="protocol",
        layer="raw_data",
    )  # Using batch_key to account for protocol differences - want to keep genes that are variable across both protocols
    sc.pp.pca(
        adata, n_comps=None, mask_var="highly_variable", random_state=0
    )  #  Roshan's workshop calculated 100 comps, then moved forward with 30, but this took long - if set to None, "Defaults to 50, or (1 - minimum dimension size of selected representation)"
    sc.pp.neighbors(
        adata, n_neighbors=30, use_rep="X_pca", metric="euclidean", random_state=0
    )  # From Roshan's workshop

    # Compute UMAP and clustering
    sc.tl.umap(adata, min_dist=0.1, random_state=0)  # From Roshan's workshop
    sc.tl.leiden(
        adata, resolution=1, random_state=0, flavor="igraph", n_iterations=2
    )  # From Roshan's workshop - scanpy recommends flavor = "igraph" and n_iterations = 2
    return