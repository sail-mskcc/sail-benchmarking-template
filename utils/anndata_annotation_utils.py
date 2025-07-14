from typing import List
import scanpy as sc

##### ANNDATA ANNOTATION UTILITIES #####
def compute_gene_category_qc_metrics(
    adatas: dict,
    category_label: str,
    pattern: str | List[str] | tuple[str],
) -> None:
    """
    Annotate genes based on a pattern and compute QC metrics for each AnnData in the dict.

    Args:
        adatas: dict with keys like 'sf_t', 'sl_t', etc. and values as AnnData objects
        category_label: name of the category to annotate (e.g., "mitochondrial", "ribosomal")
        pattern: prefix string, tuple of prefixes, or list of prefixes to match gene names
    """
    for adata in adatas.values():
        adata.var[category_label] = adata.var_names.str.startswith(pattern)
        sc.pp.calculate_qc_metrics(
            adata,
            qc_vars=[category_label],
            percent_top=None,
            log1p=False,
            inplace=True,
        )
    return


def add_gene_set_content_to_adatas(
    adatas: dict, gene_set: list, gene_set_name: str
) -> None:
    """
    Annotate and compute percent content for a specific gene set (e.g. apoptosis or housekeeping)
    across all AnnData objects in a dictionary.

    Args:
        adatas: dict of AnnData objects keyed by sample label
        gene_set: list of gene symbols representing the gene set
        gene_set_name: string label for the gene set (e.g., 'apoptosis', 'housekeeping')
    """
    for adata in adatas.values():
        filtered_genes = list(set(gene_set) & set(adata.var_names))
        print(
            f"{adata.obs.shape[0]} cells | {len(filtered_genes)} {gene_set_name} genes found"
        )

        # Mark genes belonging to the gene set
        adata.var[gene_set_name] = adata.var_names.isin(filtered_genes)

        # Calculate percent content using scanpy
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=[gene_set_name], percent_top=None, log1p=False, inplace=True
        )

    return