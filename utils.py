from typing import Callable, Dict, List, Union
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc

##### ADATA UTILITIES #####

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


##### PLOTTING UTILITIES #####

### BAR PLOTS ###

# Helper Functions
def adaptive_formatter(x, _):
    """
    Formats large numeric values for axis tick labels using adaptive human-readable notation.

    - Values ≥ 1,000,000 are shown with 'M' (e.g., 2,500,000 → '2.5M')
    - Values < 1,000,000 are shown as plain integers (e.g., 42,000 → '42000')
    """
    if x >= 1e6:
        return f"{x / 1e6:.1f}M"
    else:
        return f"{int(x)}"


# AnnData Plotting Functions
def plot_scalar_metric(
    data_dict: dict,
    metadata: dict,
    metric_func: Callable,
    metric_label: str,
    protocol_color_palette: dict,
    title: str = None,
):
    """
    Generic function to plot scalar metrics computed on each sample's data (AnnData or DataFrame),
    grouped by tissue and protocol.

    Args:
        data_dict: dict of sample data (AnnData or DataFrames), keyed by sample ID
        metadata: dict mapping sample ID to (tissue, protocol)
        metric_func: function that computes scalar metric for each sample's data
        metric_label: string for Y-axis label and plot title
        protocol_color_palette: dict mapping protocol to colors
    """
    # Assemble data
    rows = []
    for key, sample in data_dict.items():
        tissue, protocol = metadata[key]
        value = metric_func(sample)
        rows.append({"Tissue": tissue, "Protocol": protocol, "Value": value})
    df = pd.DataFrame(rows)

    # Setup
    unique_tissues = df["Tissue"].unique()
    sns.set_theme(style="white", font_scale=1.1)
    fig, axes = plt.subplots(1, len(unique_tissues), figsize=(5 * len(unique_tissues), 5), sharey=True)
    if len(unique_tissues) == 1:
        axes = [axes]

    # Plot
    for i, tissue in enumerate(unique_tissues):
        ax = axes[i]
        sns.barplot(
            data=df[df["Tissue"] == tissue],
            x="Protocol",
            y="Value",
            hue="Protocol",
            palette=protocol_color_palette,
            legend=False,
            ax=ax,
        )
        ax.set_title(tissue, fontsize=13, weight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(metric_label if i == 0 else "")
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(adaptive_formatter))
        for p in ax.patches:
            val = p.get_height()
            ax.annotate(
                f"{int(val):,}",
                (p.get_x() + p.get_width() / 2.0, val),
                ha="center", va="bottom", fontsize=10, weight="bold"
            )
    if not title:
        title = f"{metric_label} by Protocol Across Tissues"
    plt.suptitle(title, fontsize=16, weight="bold")

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    legend_elements = [
        mpl.patches.Patch(facecolor=color, label=label)
        for label, color in protocol_color_palette.items()
    ]
    fig.legend(
        handles=legend_elements,
        title="Protocol",
        loc="upper right",
        bbox_to_anchor=(1.05, 1.05),
    )
    return fig


# Function to plot a metric distribution from multiple AnnData objects and create a subplot for each unique tissue.


def plot_adata_metric_histogram(
    adatas: Dict[str, sc.AnnData],
    adata_metadata: Dict[str, tuple],
    field: str,
    axis: str = "obs",
    bins: int = 100,
    x_label: str = None,
    log_x: bool = True,
    log_y: bool = False,
    title: str = "Distribution by Protocol Across Tissues",
    row_label: str = "Cell",
    protocol_color_palette: Dict[str, str] = None,
    proportion: bool = False,
) -> None:
    """
    Plots the distribution of a field (column) from each AnnData object, grouped by tissue and protocol.
    """
    # Gather data
    rows = []
    for key, adata in adatas.items():
        tissue, protocol = adata_metadata[key]
        values = getattr(adata, axis)[field]
        if log_x:
            values = np.log1p(values)
        for val in values:
            rows.append({"Tissue": tissue, "Protocol": protocol, "Value": val})
    df = pd.DataFrame(rows)

    # Unique tissues
    unique_tissues = df["Tissue"].unique()

    # Plot
    sns.set_theme(style="white", font_scale=1.1)
    fig, axes = plt.subplots(
        1, len(unique_tissues), figsize=(5 * len(unique_tissues), 5), sharey=True
    )
    if len(unique_tissues) == 1:
        axes = [axes]

    stat_type = "probability" if proportion else "count"

    for i, tissue in enumerate(unique_tissues):
        ax = axes[i]
        sns.histplot(
            data=df[df["Tissue"] == tissue],
            x="Value",
            hue="Protocol",
            palette=protocol_color_palette,
            bins=bins,
            element="step",
            stat=stat_type,
            common_norm=False,
            legend=False,
            ax=ax,
        )
        if log_y:
            ax.set_yscale("log")

        ax.set_title(tissue, fontsize=13, weight="bold")
        ax.set_xlabel(x_label or (f"log1p({field})" if log_x else field))
        if i == 0:
            ax.set_ylabel(f"{row_label} {'Proportion' if proportion else 'Count'}")
        else:
            ax.set_ylabel("")

    # Axis scaling
    if log_x:
        raw_values = np.concatenate(
            [getattr(adatas[key], axis)[field].to_numpy() for key in adatas]
        )
        raw_values = raw_values[raw_values > 0]
        vmin, vmax = raw_values.min(), raw_values.max()
        min_exp, max_exp = int(np.floor(np.log10(vmin))), int(np.ceil(np.log10(vmax)))
        xticks_raw = [10**e for e in range(min_exp, max_exp + 1)]
        xticks_log = np.log1p(xticks_raw)

        for ax in axes:
            ax.set_xlim(np.log1p(vmin), np.log1p(vmax))
            ax.set_xticks(xticks_log)
            ax.set_xticklabels([f"{x:,}" for x in xticks_raw])
    else:
        all_values = np.concatenate(
            [getattr(adatas[key], axis)[field].to_numpy() for key in adatas]
        )
        vmin, vmax = all_values.min(), all_values.max()
        for ax in axes:
            ax.set_xlim(vmin, vmax)

    # Finalize
    plt.suptitle(title, fontsize=16, weight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    legend_elements = [
        mpl.patches.Patch(facecolor=color, label=label)
        for label, color in protocol_color_palette.items()
    ]
    fig.legend(
        handles=legend_elements,
        title="Protocol",
        loc="upper right",
        bbox_to_anchor=(1.05, 1.05),
    )
    return fig

# Function to plot a metric distribution from multiple AnnData objects and create a subplot for each unique tissue. Same as histogram but uses
# violin plots.


def plot_adata_metric_violin(
    adatas: Dict[str, sc.AnnData],
    adata_metadata: Dict[str, tuple],
    field: str,
    axis: str = "obs",
    log_y: bool = True,
    title: str = "Distribution by Protocol Across Tissues",
    y_label: str = "Value",
    protocol_color_palette: Dict[str, str] = None,
) -> None:
    """
    Plots a violin distribution of a given .obs or .var field from each AnnData object,
    grouped by tissue and protocol. One subplot is generated per tissue.

    Args:
        adatas: Dictionary of AnnData objects keyed by sample ID (e.g., 'sf_ln').
        adata_metadata: Dictionary mapping sample ID to (tissue, protocol) tuples.
        field: Column name to plot from .obs or .var (e.g., 'total_counts').
        axis: Whether to pull from .obs or .var (default: 'obs').
        log_y: Whether to apply log1p transformation to the y-axis values (default: True).
        title: Title of the entire figure.
        y_label: Label for the y-axis.
        protocol_color_palette: Optional dictionary mapping protocol names to colors.
    """
    # Assemble data
    rows = []
    for key, adata in adatas.items():
        tissue, protocol = adata_metadata[key]
        values = getattr(adata, axis)[field]
        if log_y:
            values = np.log1p(values)
        for val in values:
            rows.append({"Tissue": tissue, "Protocol": protocol, "Value": val})
    df = pd.DataFrame(rows)

    # Get unique tissues
    unique_tissues = df["Tissue"].unique()

    # Setup plot grid
    sns.set_theme(style="white", font_scale=1.1)
    fig, axes = plt.subplots(
        1, len(unique_tissues), figsize=(5 * len(unique_tissues), 5), sharey=True
    )
    if len(unique_tissues) == 1:
        axes = [axes]

    for i, tissue in enumerate(unique_tissues):
        ax = axes[i]
        sns.violinplot(
            data=df[df["Tissue"] == tissue],
            x="Protocol",
            y="Value",
            hue="Protocol",
            palette=protocol_color_palette,
            ax=ax,
            linewidth=1,
        )
        ax.set_title(tissue, fontsize=13, weight="bold")
        ax.set_xlabel("")
        ax.set_ylabel(y_label if i == 0 else "")

    plt.suptitle(title, fontsize=16, weight="bold")

    # Custom y-axis ticks if log-transformed
    if log_y:
        raw_values = np.concatenate(
            [getattr(adatas[key], axis)[field].to_numpy() for key in adatas]
        )
        raw_values = raw_values[raw_values > 0]

        ymin = raw_values.min()
        ymax = raw_values.max()
        min_exp = int(np.floor(np.log10(ymin)))
        max_exp = int(np.ceil(np.log10(ymax)))
        yticks_raw = [10**e for e in range(min_exp, max_exp + 1)]
        yticks_log = np.log1p(yticks_raw)

        for ax in axes:
            ax.set_ylim(np.log1p(ymin), np.log1p(ymax))
            ax.set_yticks(yticks_log)
            ax.set_yticklabels([f"{y:,}" for y in yticks_raw])

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


def plot_celltype_proportions_by_protocol(
    adata: sc.AnnData,
    tissue: str,
    protocol_color_palette: dict = None,
    annotate: bool = False,
    ) -> plt.Figure:
    """
    Plots bar chart of cell type proportions by protocol for a given tissue using Seaborn.

    Args:
        adata: AnnData object for the specific tissue.
        tissue: Tissue name to extract and plot.
        annotate: Whether to annotate bars with raw counts.
        protocol_color_palette: Optional dictionary mapping protocol to color.

    Returns:
        Matplotlib figure object.
    """
    obs = adata.obs[["protocol", "majority_voting"]].copy()

    # Count and proportion
    counts = (
        obs.groupby(["protocol", "majority_voting"]).size().reset_index(name="count")
    )
    total_per_protocol = counts.groupby("protocol")["count"].transform("sum")
    counts["proportion"] = counts["count"] / total_per_protocol

    # Sort cell types
    counts["majority_voting"] = counts["majority_voting"].astype(str)
    counts = counts.sort_values("majority_voting")
    counts["percentage"] = counts["proportion"] * 100

    # Set up the plot
    sns.set_theme(style="white", font_scale=1.1)
    fig, ax = plt.subplots(figsize=(10, 6))

    sns.barplot(
        data=counts,
        x="majority_voting",
        y="percentage",
        hue="protocol",
        palette=protocol_color_palette,
        ax=ax,
        dodge=True,
    )

    if annotate:
        for p in ax.patches:
            height = p.get_height()
            if height > 0:
                ax.annotate(
                    f"{height:.2f}",
                    (p.get_x() + p.get_width() / 2.0, height),
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )

    # Format axes and legend
    ax.set_title(
        f"Cell Type Proportions by Protocol ({tissue})", fontsize=14, weight="bold"
    )
    ax.set_xlabel("Cell Type", fontsize=12)
    ax.set_ylabel("Percentage", fontsize=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.legend(title="Protocol", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    return fig
    

def plot_cluster_protocol_stackplots(
    combined_by_tissue: Dict[str, sc.AnnData],
    cluster_key: str = "leiden",
    protocol_key: str = "protocol",
    protocol_color_palette: Dict[str, str] = None,
    tissue_order: List[str] = None,
    title: str = None,
) -> plt.Figure:
    """
    For each tissue, plot a stacked barplot showing how many cells from each protocol
    make up each cluster (cluster_key).

    Args:
        combined_by_tissue (dict): Dict of AnnData objects per tissue name.
        cluster_key (str): Column in .obs for clustering (e.g. 'leiden').
        protocol_key (str): Column in .obs for protocol (e.g. 'protocol').
        protocol_palette (dict): Mapping of protocol to color.
        tissue_order (list): Optional order of tissues to plot.

    Returns:
        fig (matplotlib.Figure): The figure containing the stacked barplots.
    """
    sns.set_theme(style="white", font_scale=1.1)
    # Build long-form dataframe with one row per (Tissue, Cluster, Protocol)
    records = []
    for tissue, adata in combined_by_tissue.items():
        df = adata.obs[[cluster_key, protocol_key]].copy()
        df["Tissue"] = tissue
        records.append(df)

    df_all = pd.concat(records)
    df_all[cluster_key] = df_all[cluster_key].astype(int)

    # Count cells per group
    counts = (
        df_all.groupby(["Tissue", cluster_key, protocol_key], observed=True)
        .size()
        .reset_index(name="Count")
    )

    # Pivot to wide form then normalize row-wise
    pivoted = counts.pivot_table(
        index=["Tissue", cluster_key],
        columns=protocol_key,
        values="Count",
        fill_value=0,
    )
    pivoted_norm = pivoted.div(pivoted.sum(axis=1), axis=0).reset_index()

    unique_tissues = tissue_order or pivoted_norm["Tissue"].unique()

    fig, axes = plt.subplots(
        1, len(unique_tissues), figsize=(6 * len(unique_tissues), 5), sharey=True
    )

    if len(unique_tissues) == 1:
        axes = [axes]

    for i, tissue in enumerate(unique_tissues):
        ax = axes[i]
        df_tissue = pivoted_norm[pivoted_norm["Tissue"] == tissue]
        clusters = df_tissue[cluster_key]

        bottom = pd.Series([0.0] * len(clusters), index=clusters)
        for protocol in protocol_color_palette.keys():
            if protocol in df_tissue.columns:
                values = df_tissue[protocol].values
                ax.bar(
                    clusters,
                    values,
                    bottom=bottom.loc[clusters].values,
                    color=protocol_color_palette[protocol],
                    label=protocol,
                )
                bottom.loc[clusters] += values

        ax.set_title(tissue, fontsize=13, weight="bold")
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Proportion of Cells" if i == 0 else "")
        ax.set_xticks(range(len(clusters)))
        ax.set_xticklabels(clusters, rotation=45)
        ax.set_ylim(0, 1)

    fig.legend(
        handles=[
            mpl.patches.Patch(color=c, label=p) for p, c in protocol_color_palette.items()
        ],
        title="Protocol",
        loc="upper right",
        bbox_to_anchor=(1.05, 1),
    )
    if not title:
        title = "Cluster Composition by Protocol Across Tissues"
    plt.suptitle(
        title, fontsize=16, weight="bold"
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig

def plot_doublet_stack_by_cluster(
    combined_by_tissue: Dict[str, "sc.AnnData"],
    tissue: str,
    color_palette: Dict[str, str],
) -> None:
    """
    Plot stacked barplot of doublet vs non-doublet counts per cluster, grouped by protocol.

    Args:
        combined_by_tissue: Dict of tissue name to AnnData.
        tissue: Tissue name to plot.
        color_palette: Dict of protocol -> hex color (opaque), e.g. {"Singulator+FACS": "#6B9EB2"}
    """
    adata = combined_by_tissue[tissue].copy()
    adata.obs["leiden"] = adata.obs["leiden"].astype(str)

    # Total cells per cluster/protocol
    cell_counts = (
        adata.obs.groupby(["leiden", "protocol"]).size().reset_index(name="total_cells")
    )

    # Doublet counts
    doublet_counts = (
        adata.obs.groupby(["leiden", "protocol"])["predicted_doublet"]
        .sum()
        .reset_index()
        .rename(columns={"predicted_doublet": "doublet_count"})
    )

    # Merge and compute non-doublets
    summary = pd.merge(doublet_counts, cell_counts, on=["leiden", "protocol"])
    summary["non_doublets"] = summary["total_cells"] - summary["doublet_count"]
    summary["leiden"] = summary["leiden"].astype(int)
    summary = summary.sort_values("leiden")

    # Plot setup
    fig = plt.figure(figsize=(9, 6))
    bar_width = 0.35
    leiden_clusters = summary["leiden"].unique()
    x = np.arange(len(leiden_clusters))
    protocols = summary["protocol"].unique()
    offsets = {protocol: (i - 0.5) * bar_width for i, protocol in enumerate(protocols)}
    faded_colors = {protocol: color_palette[protocol] + "80" for protocol in protocols}

    for protocol in protocols:
        subset = summary[summary["protocol"] == protocol]
        indices = [np.where(leiden_clusters == cl)[0][0] for cl in subset["leiden"]]

        plt.bar(
            x[indices] + offsets[protocol],
            subset["non_doublets"],
            width=bar_width,
            color=faded_colors[protocol],
            label=f"{protocol} (non-doublets)",
        )
        plt.bar(
            x[indices] + offsets[protocol],
            subset["doublet_count"],
            width=bar_width,
            bottom=subset["non_doublets"],
            color=color_palette[protocol],
            label=f"{protocol} (doublets)",
        )

    plt.xticks(x, leiden_clusters, rotation=45)
    plt.xlabel("Cluster")
    plt.ylabel("Cell Count")
    plt.title(f"Doublets vs Non-Doublets per Cluster by Protocol ({tissue})")
    plt.legend(loc="upper right", title="Protocol + Cell Type")
    plt.tight_layout()
    return fig

def plot_cluster_metric_boxplots(
    combined_by_tissue: Dict[str, sc.AnnData],
    metric: str = "total_counts",
    cluster_key: str = "leiden",
    tissue_order: List[str] = None,
    title: str = None,
    log_scale: bool = False,
    axline: int = None,
    vertical: bool = False,
) -> plt.Figure:
    """
    For each tissue, plot boxplots of a given QC metric per cluster (e.g., total_counts).

    Args:
        combined_by_tissue (dict): Dict of AnnData objects per tissue name.
        metric (str): The .obs column to plot (e.g., 'total_counts').
        cluster_key (str): Column in .obs for clustering (e.g. 'leiden').
        tissue_order (list): Optional order of tissues to plot.
        title (str): Optional overall plot title.
        log_scale (bool): Whether to use log scale for the metric axis.
        axline (int): Optional reference line.
        vertical (bool): If True, boxplots are horizontal (metric on x-axis).
    """
    sns.set_theme(style="white", font_scale=1.1)

    # Build a long-form DataFrame: Tissue, Cluster, Metric
    records = []
    for tissue, adata in combined_by_tissue.items():
        df = adata.obs[[cluster_key, metric]].copy()
        df["Tissue"] = tissue
   
        df[cluster_key] = df[cluster_key].astype(int)
        records.append(df)

    df_all = pd.concat(records)
    unique_tissues = tissue_order or df_all["Tissue"].unique()

    # Plot setup
    if vertical:
        fig, axes = plt.subplots(
            len(unique_tissues), 1, figsize=(7, 6 * len(unique_tissues)), sharex=True
        )
    else:
        fig, axes = plt.subplots(
            1, len(unique_tissues), figsize=(6 * len(unique_tissues), 5), sharey=True
        )

    if len(unique_tissues) == 1:
        axes = [axes]

    for i, tissue in enumerate(unique_tissues):
        ax = axes[i]
        df_tissue = df_all[df_all["Tissue"] == tissue]

        if vertical:
            # Convert cluster to string for categorical handling
            df_tissue = df_tissue.copy()
            df_tissue[cluster_key] = df_tissue[cluster_key].astype(str)
            
            # Get cluster order as strings
            cluster_order = sorted(
                df_tissue[cluster_key].unique(),
                key=lambda x: int(x) if x.isdigit() else x
            )
            
            sns.boxplot(
                data=df_tissue,
                x=metric,
                y=cluster_key,
                ax=ax,
                color="#5fbced",
                fliersize=2,
                linewidth=0.75,
                order=cluster_order,
            )

            ax.margins(y=0.01)

            # Apply log scale to x-axis if requested
            if log_scale:
                ax.set_xscale('log')
                ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(adaptive_formatter))
            
            ax.set_xlabel(metric.replace("_", " ").title())
            ax.set_ylabel("Cluster")
            
            # Add vertical reference line if specified
            if axline is not None:
                ax.axvline(axline, color="red", linestyle="--", linewidth=1)
        else:
            # Convert cluster to string for categorical handling
            df_tissue = df_tissue.copy()
            df_tissue[cluster_key] = df_tissue[cluster_key].astype(str)
            
            # Get cluster order as strings
            cluster_order = sorted(
                df_tissue[cluster_key].unique(),
                key=lambda x: int(x) if x.isdigit() else x
            )
            
            sns.boxplot(
                data=df_tissue,
                x=cluster_key,
                y=metric,
                ax=ax,
                color="#5fbced",
                fliersize=2,
                linewidth=0.75,
                order=cluster_order,
            )
            
            # Apply log scale to y-axis if requested
            if log_scale:
                ax.set_yscale('log')
                ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(adaptive_formatter))

            ax.set_xlabel("Cluster")
            ax.set_ylabel(metric.replace("_", " ").title())
            ax.tick_params(axis='x', rotation=45)
            
            # Add horizontal reference line if specified
            if axline is not None:
                ax.axhline(axline, color="red", linestyle="--", linewidth=1)

        ax.set_title(tissue, fontsize=13, weight="bold")

    if not title:
        title = f"{metric.replace('_', ' ').title()} Distribution by Cluster Across Tissues"
    plt.suptitle(title, fontsize=16, weight="bold")
    plt.tight_layout(rect=[0, 0, 0.95, 0.99])
    return fig


### UMAP PLOTS ###

def plot_umap_by_obs_feature(
    combined_by_tissue: Dict[str, sc.AnnData],
    feature: str,
    color_palette: Union[Dict[str, str], str] = "viridis",
    title: str = None,
    size: float = 8.0,
    alpha: float = 0.8,
    log_scale: bool = False,
    clip_values: tuple = None,
    shuffle: bool = True,
    vertical: bool = False,
    legend_bbox_to_anchor: tuple = (1.05, 1.05)
) -> plt.Figure:
    """
    Plot UMAPs for each tissue colored by a feature in .obs.

    Args:
        combined_by_tissue: Dictionary of AnnData objects per tissue.
        feature: Column in .obs to plot (categorical or continuous).
        color_palette: Dict for categorical or colormap name for continuous.
        title: Figure title.
        size: Dot size.
        alpha: Dot transparency.
        log_scale: Use log1p scale (applies to continuous features).
        clip_values: (low, high) quantiles to clip (e.g., (0.01, 0.99)).
        shuffle: Shuffle points before plotting (default True).

    Returns:
        Matplotlib figure.
    """
    sns.set_theme(style="white")
    tissues = list(combined_by_tissue.keys())
    if vertical:
        fig, axes = plt.subplots(
            len(tissues), 1, figsize=(6, 5 * len(tissues)), squeeze=False
        )
    else:
        fig, axes = plt.subplots(
            1, len(tissues), figsize=(6 * len(tissues), 5), squeeze=False
        )

    # Determine global vmin/vmax
    all_values = pd.concat(
        [adata.obs[feature] for adata in combined_by_tissue.values()]
    )
    clip_low = (
        max(1, all_values.quantile(clip_values[0]))
        if clip_values and log_scale
        else all_values.min()
    )
    clip_high = all_values.quantile(clip_values[1]) if clip_values else all_values.max()
    vmin, vmax = (
        (np.log1p(clip_low), np.log1p(clip_high))
        if log_scale
        else (clip_low, clip_high)
    )

    # Convert string colormap to dict for categorical legend if needed
    def make_categorical_palette(vals, cmap_name):
        cmap = mpl.colormaps.get_cmap(cmap_name).resampled(len(vals))
        return {val: mpl.colors.to_hex(cmap(i)) for i, val in enumerate(sorted(vals))}

    for i, tissue in enumerate(tissues):
        ax = axes[i, 0] if vertical else axes[0, i]
        adata = combined_by_tissue[tissue]
        coords = adata.obsm["X_umap"]
        values = adata.obs[feature].copy()
        # Sort to plot False/0 first, True/1 last (on top)
        if pd.api.types.is_bool_dtype(values) or set(values.unique()) <= {0, 1}:
            values = values.astype(bool)
            sorted_idx = np.argsort(values.values)  # False first, True last
            coords = coords[sorted_idx]
            values = values.iloc[sorted_idx]
        elif shuffle:
            np.random.seed(0)
            shuffled = np.random.permutation(adata.n_obs)
            coords = coords[shuffled]
            values = values.iloc[shuffled]

        if (
            values.dtype.name == "category"
            or values.dtype == object
            or isinstance(color_palette, dict)
        ):
            # Categorical
            if isinstance(color_palette, str):
                unique_vals = sorted(values.dropna().unique())
                color_palette = make_categorical_palette(unique_vals, color_palette)
            sns.scatterplot(
                x=coords[:, 0],
                y=coords[:, 1],
                hue=values,
                palette=color_palette,
                ax=ax,
                s=size,
                alpha=alpha,
                linewidth=0,
            )
            ax.legend_.remove()
        else:
            # Continuous
            values = values.clip(clip_low, clip_high)
            if log_scale:
                values = np.log1p(values)
            scatter = ax.scatter(
                coords[:, 0],
                coords[:, 1],
                c=values,
                cmap=color_palette,
                vmin=vmin,
                vmax=vmax,
                s=size,
                alpha=alpha,
                linewidths=0,
            )
            if i == len(tissues) - 1:
                if vertical:
                    cb = plt.colorbar(scatter, ax=ax, location = "bottom", orientation = "horizontal", fraction=0.046, pad=0.15)
                else:
                    cb = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.14)
                if log_scale:
                    raw_ticks = np.logspace(
                        np.log10(max(1, clip_low)), np.log10(clip_high), num=6
                    )
                    cb.set_ticks(np.log1p(raw_ticks))
                    cb.set_ticklabels(
                        [
                            f"{int(t/1000)}k" if t >= 1000 else str(int(t))
                            for t in raw_ticks
                        ]
                    )
                else:
                    raw_ticks = np.linspace(clip_low, clip_high, 6)
                    cb.set_ticks(raw_ticks)
                    cb.set_ticklabels([f"{int(t):,}" for t in raw_ticks])
                
                rotation = 0 if vertical else 270
                labelpad = -160 if vertical else 15
                loc = "top" if vertical else "center"

                cb.ax.set_ylabel(
                    feature.replace("log1p_", "").replace("_", " ").title(),
                    rotation=rotation,
                    labelpad=labelpad,
                    loc=loc,
                )

        ax.set_title(tissue, weight="bold")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2" if i == 0 else "")
        ax.set_xticks([])
        ax.set_yticks([])

    # Add global legend for categorical data
    if (
        values.dtype.name == "category"
        or values.dtype == object
        or pd.api.types.is_bool_dtype(values)
    ):
        unique_vals = sorted(
            values.dropna().unique(),
            key=lambda x: float(x) if str(x).replace(".", "", 1).isdigit() else str(x),
        )
        legend_elements = [
            mpl.patches.Patch(facecolor=color_palette.get(val, "gray"), label=val)
            for val in unique_vals
        ]
        if vertical:
            fig.legend(
                handles=legend_elements,
                title=feature,
                loc="upper right",
                bbox_to_anchor=legend_bbox_to_anchor,
                fontsize="small",
                title_fontsize="medium",
            )
        else:
            fig.legend(
                handles=legend_elements,
                title=feature,
                loc="upper right",
                bbox_to_anchor=legend_bbox_to_anchor,
                fontsize="small",
                title_fontsize="medium",
            )

    if not title:
        title = f"UMAP Colored by {feature.replace('_', ' ').title()}"
        if log_scale:
            title += " (Log Scale)"
        if clip_values:
            title += f" (Clipped {clip_values[0]*100:.1f}%-{clip_values[1]*100:.1f}%)"

    plt.suptitle(title, fontsize=16, weight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig