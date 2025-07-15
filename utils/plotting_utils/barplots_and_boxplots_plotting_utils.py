from typing import Dict, List, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from .common_plotting_utils import adaptive_formatter


def plot_celltype_proportions_by_protocol(
    adata: sc.AnnData,
    tissue: str,
    annotate: bool = False,
    metric: str = "majority_voting",
    protocol_key: str = "protocol",
    protocol_color_palette: dict = None,
    percentage: bool = False,
    plot_width: float = 10.0,
    plot_height: float = 6.0,
    legend_bbox_to_anchor: tuple = (1.05, 1),
    title: str = None,
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
    obs = adata.obs[[protocol_key, metric]].copy()

    # Count and proportion
    counts = obs.groupby([protocol_key, metric]).size().reset_index(name="count")
    total_per_protocol = counts.groupby(protocol_key)["count"].transform("sum")
    counts["proportion"] = counts["count"] / total_per_protocol
    counts["percentage"] = counts["proportion"] * 100

    # Sort cell types
    counts[metric] = counts[metric].astype(str)
    counts = counts.sort_values(metric)

    # Set up the plot
    sns.set_theme(style="white", font_scale=1.1)
    fig, ax = plt.subplots(figsize=(plot_width, plot_height))

    # Use the provided color palette or default to a seaborn palette
    if protocol_color_palette is None:
        protocol_color_palette = sns.color_palette(
            "husl", n_colors=counts[protocol_key].nunique()
        )
        protocol_color_palette = dict(
            zip(counts[protocol_key].unique(), protocol_color_palette)
        )

    # Plotting
    value_to_plot = "percentage" if percentage else "count"
    sns.barplot(
        data=counts,
        x=metric,
        y=value_to_plot,
        hue=protocol_key,
        palette=protocol_color_palette,
        ax=ax,
        dodge=True,
    )

    # Optionally annotate bars
    if annotate:
        for p in ax.patches:
            height = p.get_height()
            if height > 0:
                annotation = f"{int(height)}" if not percentage else f"{height:.2f}%"
                ax.annotate(
                    annotation,
                    (p.get_x() + p.get_width() / 2.0, height),
                    ha="center",
                    va="bottom",
                    fontsize=9,
                )

    # Set title and labels
    if title is None:
        title = (
            f"Cell Type Percentages by Protocol ({tissue})"
            if percentage
            else f"Cell Type Counts by Protocol ({tissue})"
        )

    ax.set_title(title, fontsize=14, weight="bold")

    # Format axes and legend
    ylabel = "Percentage" if percentage else "Cell Type Count"
    ax.set_xlabel("Cell Type", fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.legend(title="Protocol", bbox_to_anchor=legend_bbox_to_anchor, loc="upper left")
    plt.tight_layout()
    return fig


def plot_doublet_stack_by_cluster(
    combined_by_tissue: Dict[str, sc.AnnData],
    tissue: str,
    cluster_key: str = "leiden",
    protocol_key: str = "protocol",
    doublet_key: str = "predicted_doublet",
    bar_width: float = 0.35,
    figsize: tuple = (9, 6),
    protocol_color_palette: Dict[str, str] = None,
) -> None:
    """
    Plot stacked barplot of doublet vs non-doublet counts per cluster, grouped by protocol.

    Args:
        combined_by_tissue: Dict of tissue name to AnnData.
        tissue: Tissue name to plot.
        color_palette: Dict of protocol -> hex color (opaque), e.g. {"Singulator+FACS": "#6B9EB2"}
    """
    adata = combined_by_tissue[tissue].copy()
    adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)

    # Total cells per cluster/protocol
    cell_counts = (
        adata.obs.groupby([cluster_key, protocol_key])
        .size()
        .reset_index(name="total_cells")
    )

    # Doublet counts
    doublet_counts = (
        adata.obs.groupby([cluster_key, protocol_key])[doublet_key]
        .sum()
        .reset_index()
        .rename(columns={doublet_key: "doublet_count"})
    )

    # Merge and compute non-doublets
    summary = pd.merge(doublet_counts, cell_counts, on=[cluster_key, protocol_key])
    summary["non_doublets"] = summary["total_cells"] - summary["doublet_count"]
    summary[cluster_key] = summary[cluster_key].astype(int)
    summary = summary.sort_values(cluster_key)

    # Plot setup
    fig = plt.figure(figsize=figsize)

    # Arrange Clusters
    leiden_clusters = summary[cluster_key].unique()
    x = np.arange(len(leiden_clusters))

    # Create offsets for each protocol
    unique_protocols = summary[protocol_key].unique()
    offsets = {
        protocol: (i - 0.5) * bar_width for i, protocol in enumerate(unique_protocols)
    }

    # Create protocol color palette if not provided
    if protocol_color_palette is None:
        protocol_color_palette = sns.color_palette(
            "husl", n_colors=len(unique_protocols)
        )
        protocol_color_palette = dict(zip(protocols, protocol_color_palette))

    faded_colors = {
        protocol: protocol_color_palette[protocol] + "80"
        for protocol in unique_protocols
    }

    for protocol in unique_protocols:
        subset = summary[summary[protocol_key] == protocol]
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
            color=protocol_color_palette[protocol],
            label=f"{protocol} (doublets)",
        )

    plt.xticks(x, leiden_clusters, rotation=45)
    plt.xlabel("Cluster")
    plt.ylabel("Cell Count")
    plt.title(f"Doublets vs Non-Doublets per Cluster by Protocol ({tissue})")
    plt.legend(loc="upper right", title="Protocol + Cell Type")
    plt.tight_layout()
    return fig


# def plot_categorical_stack_by_cluster(
#     combined_by_tissue: Dict[str, sc.AnnData],
#     tissue: str,
#     category_key: str,
#     cluster_key: str = "leiden",
#     protocol_key: str = "protocol",
#     bar_width: float = 0.35,
#     figsize: tuple = (9, 6),
#     protocol_color_palette: Dict[str, str] = None,
#     category_color_palette: Dict[str, str] = None,
# ) -> plt.Figure:
#     """
#     Plot stacked barplot of categorical counts per cluster, grouped by protocol.

#     Args:
#         combined_by_tissue: Dict of tissue name to AnnData.
#         tissue: Tissue name to plot.
#         category_key: Key in `.obs` for the categorical variable (e.g., "predicted_doublet").
#         cluster_key: Cluster identifier key in `.obs` (e.g., "leiden").
#         protocol_key: Protocol identifier key in `.obs`.
#         bar_width: Width of each bar.
#         figsize: Tuple for figure size.
#         protocol_color_palette: Dict of protocol -> color.
#         category_color_palette: Dict of category -> color.
#     """
#     import matplotlib.pyplot as plt
#     import seaborn as sns
#     import numpy as np
#     import pandas as pd

#     adata = combined_by_tissue[tissue].copy()
#     adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)

#     # Get counts per category/cluster/protocol
#     summary = (
#         adata.obs
#         .groupby([cluster_key, protocol_key, category_key])
#         .size()
#         .reset_index(name="count")
#     )

#     # Ensure consistent cluster ordering
#     summary[cluster_key] = summary[cluster_key].astype(int)
#     summary = summary.sort_values(cluster_key)
#     clusters = summary[cluster_key].unique()
#     x = np.arange(len(clusters))

#     # Get unique protocols and categories
#     protocols = summary[protocol_key].unique()
#     categories = summary[category_key].unique()

#     # Bar offset per protocol
#     offsets = {p: (i - 0.5) * bar_width for i, p in enumerate(protocols)}

#     # Protocol color palette
#     if protocol_color_palette is None:
#         protocol_color_palette = dict(zip(
#             protocols, sns.color_palette("Set2", len(protocols))
#         ))

#     # Category color palette (within each bar)
#     if category_color_palette is None:
#         category_color_palette = dict(zip(
#             categories, sns.color_palette("pastel", len(categories))
#         ))

#     # Plot
#     fig = plt.figure(figsize=figsize)

#     for protocol in protocols:
#         proto_data = summary[summary[protocol_key] == protocol]
#         grouped = proto_data.groupby(cluster_key)

#         for cluster in clusters:
#             cluster_data = grouped.get_group(cluster) if cluster in grouped.groups else pd.DataFrame()
#             bottom = 0
#             for category in categories:
#                 value = cluster_data[cluster_data[category_key] == category]["count"]
#                 count = int(value.values[0]) if not value.empty else 0
#                 plt.bar(
#                     x[clusters == cluster] + offsets[protocol],
#                     [count],
#                     width=bar_width,
#                     bottom=bottom,
#                     color=category_color_palette[category],
#                     label=f"{protocol} ({category})" if cluster == clusters[0] else None,
#                 )
#                 bottom += count

#     plt.xticks(x, clusters, rotation=45)
#     plt.xlabel("Cluster")
#     plt.ylabel("Cell Count")
#     plt.title(f"{category_key} Distribution per Cluster by Protocol ({tissue})")
#     handles, labels = plt.gca().get_legend_handles_labels()
#     by_label = dict(zip(labels, handles))
#     plt.legend(by_label.values(), by_label.keys(), title="Protocol + Category", loc="upper right")
#     plt.tight_layout()
#     return fig


def plot_categorical_stack_by_cluster(
    combined_by_tissue: Dict[str, sc.AnnData],
    tissue: str,
    category_key: str,
    cluster_key: str = "leiden",
    protocol_key: str = "protocol",
    bar_width: float = 0.35,
    figsize: tuple = (9, 6),
    protocol_color_palette: Dict[str, str] = None,
    highlight_false: bool = False,
) -> plt.Figure:
    """
    Plot stacked barplot of categorical counts per cluster, grouped by protocol.

    Args:
        combined_by_tissue: Dict of tissue name to AnnData.
        tissue: Tissue name to plot.
        category_key: Key in `.obs` for the categorical variable (e.g., "CellBender_Included").
        cluster_key: Cluster identifier key in `.obs` (e.g., "leiden").
        protocol_key: Protocol identifier key in `.obs`.
        bar_width: Width of each bar.
        figsize: Tuple for figure size.
        protocol_color_palette: Dict of protocol -> base hex color.
    """

    adata = combined_by_tissue[tissue].copy()
    adata.obs[cluster_key] = adata.obs[cluster_key].astype(str)

    # Get counts per category/cluster/protocol
    summary = (
        adata.obs.groupby([cluster_key, protocol_key, category_key])
        .size()
        .reset_index(name="count")
    )

    summary[cluster_key] = summary[cluster_key].astype(int)
    summary = summary.sort_values(cluster_key)
    clusters = summary[cluster_key].unique()
    x = np.arange(len(clusters))

    # Get unique protocols and categories
    protocols = summary[protocol_key].unique()
    categories = summary[category_key].unique()

    # Assign offsets
    offsets = {p: (i - 0.5) * bar_width for i, p in enumerate(protocols)}

    # Protocol base colors
    if protocol_color_palette is None:
        protocol_color_palette = dict(
            zip(protocols, sns.color_palette("Set2", len(protocols)))
        )

    alpha_values = (1.0, 0.4)
    if highlight_false:
        alpha_values = (0.4, 1.0)  # Highlight 'False' category

    protocol_category_palette = {
        (protocol, category): mpl.colors.to_rgba(
            protocol_color_palette[protocol],
            alpha=alpha_values[0] if bool(category) else alpha_values[1],
        )
        for protocol in protocols
        for category in categories
    }
    # Plot
    fig = plt.figure(figsize=figsize)

    for protocol in protocols:
        proto_data = summary[summary[protocol_key] == protocol]
        grouped = proto_data.groupby(cluster_key)

        for cluster in clusters:
            cluster_data = (
                grouped.get_group(cluster)
                if cluster in grouped.groups
                else pd.DataFrame()
            )
            bottom = 0
            for category in categories:
                value = cluster_data[cluster_data[category_key] == category]["count"]
                count = int(value.values[0]) if not value.empty else 0
                plt.bar(
                    x[clusters == cluster] + offsets[protocol],
                    [count],
                    width=bar_width,
                    bottom=bottom,
                    color=protocol_category_palette[(protocol, category)],
                    label=(
                        f"{protocol} ({category})" if cluster == clusters[0] else None
                    ),
                )
                bottom += count

    plt.xticks(x, clusters, rotation=45)
    plt.xlabel("Cluster")
    plt.ylabel("Cell Count")
    plt.title(f"{category_key} Distribution per Cluster by Protocol ({tissue})")

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(
        by_label.values(),
        by_label.keys(),
        title="Protocol + Category",
        loc="upper right",
    )
    plt.tight_layout()
    return fig


def plot_cluster_protocol_stackplots(
    combined_by_tissue: Dict[str, sc.AnnData],
    cluster_key: str = "leiden",
    protocol_key: str = "protocol",
    protocol_color_palette: Dict[str, str] = None,
    tissue_order: List[str] = None,
    tissue_key: str = "Tissue",
    title: str = None,
    subplot_width: float = 6.0,
    subplot_height: float = 5.0,
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
        df[tissue_key] = tissue
        records.append(df)

    df_all = pd.concat(records)
    df_all[cluster_key] = df_all[cluster_key].astype(int)

    # Count cells per group
    counts = (
        df_all.groupby([tissue_key, cluster_key, protocol_key], observed=True)
        .size()
        .reset_index(name="Count")
    )

    # Pivot to wide form then normalize row-wise
    pivoted = counts.pivot_table(
        index=[tissue_key, cluster_key],
        columns=protocol_key,
        values="Count",
        fill_value=0,
    )

    # Normalize to proportions
    pivoted_norm = pivoted.div(pivoted.sum(axis=1), axis=0).reset_index()

    # Get unique tissues and order them if specified
    unique_tissues = tissue_order or pivoted_norm[tissue_key].unique()

    # Get unique protocols and create color palette if not provided
    if protocol_color_palette is None:
        protocols = pivoted_norm[protocol_key].unique()
        protocol_color_palette = sns.color_palette("husl", n_colors=len(protocols))
        protocol_color_palette = dict(zip(protocols, protocol_color_palette))

    # Set up the plot
    fig, axes = plt.subplots(
        1,
        len(unique_tissues),
        figsize=(subplot_width * len(unique_tissues), subplot_height),
        sharey=True,
    )

    if len(unique_tissues) == 1:
        axes = [axes]

    # Plot each tissue's stacked barplot
    for i, tissue in enumerate(unique_tissues):
        ax = axes[i]
        df_tissue = pivoted_norm[pivoted_norm[tissue_key] == tissue]
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

    # Add legend and title
    fig.legend(
        handles=[
            mpl.patches.Patch(color=c, label=p)
            for p, c in protocol_color_palette.items()
        ],
        title="Protocol",
        loc="upper right",
        bbox_to_anchor=(1.05, 1),
    )
    if not title:
        title = "Cluster Composition by Protocol Across Tissues"
    plt.suptitle(title, fontsize=16, weight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig


def plot_cluster_metric_boxplots(
    combined_by_tissue: Dict[str, sc.AnnData],
    metric: str = "total_counts",
    cluster_key: str = "leiden",
    tissue_key: str = "Tissue",
    tissue_order: List[str] = None,
    title: str = None,
    log_scale: bool = False,
    axline: int = None,
    vertical: bool = False,
    fig_width: float = None,
    fig_height: float = None,
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
        df[tissue_key] = tissue

        # df[cluster_key] = df[cluster_key].astype(int)
        if pd.api.types.is_numeric_dtype(df[cluster_key]):
            df[cluster_key] = df[cluster_key].astype(int)
        records.append(df)

    df_all = pd.concat(records)
    if tissue_order is None:
        unique_tissues = df_all[tissue_key].unique()
    else:
        unique_tissues = list(set(tissue_order))

    # Plot setup
    if vertical:
        fig, axes = plt.subplots(
            len(unique_tissues),
            1,
            figsize=(fig_width or 7, (fig_height or 6) * len(unique_tissues)),
            sharex=True,
        )
    else:
        fig, axes = plt.subplots(
            1,
            len(unique_tissues),
            figsize=(fig_width or 6 * len(unique_tissues), fig_height or 5),
            sharey=True,
        )

    if len(unique_tissues) == 1:
        axes = [axes]

    for i, tissue in enumerate(unique_tissues):
        ax = axes[i]
        df_tissue = df_all[df_all["Tissue"] == tissue]
        df_tissue[cluster_key] = df_tissue[cluster_key].astype(str)

        # Get cluster order as strings
        cluster_order = sorted(
            df_tissue[cluster_key].unique(),
            key=lambda x: int(x) if x.isdigit() else x,
        )

        if vertical:
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
                ax.set_xscale("log")
                ax.xaxis.set_major_formatter(
                    mpl.ticker.FuncFormatter(adaptive_formatter)
                )

            ax.set_xlabel(metric.replace("_", " ").title())
            ax.set_ylabel("Cluster")

            # Add vertical reference line if specified
            if axline is not None:
                ax.axvline(axline, color="red", linestyle="--", linewidth=1)
        else:
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
                ax.set_yscale("log")
                ax.yaxis.set_major_formatter(
                    mpl.ticker.FuncFormatter(adaptive_formatter)
                )

            ax.set_xlabel("Cluster")
            ax.set_ylabel(metric.replace("_", " ").title())
            ax.tick_params(axis="x", rotation=45)

            # Add horizontal reference line if specified
            if axline is not None:
                ax.axhline(axline, color="red", linestyle="--", linewidth=1)

        ax.set_title(tissue, fontsize=13, weight="bold")

    if not title:
        title = (
            f"{metric.replace('_', ' ').title()} Distribution by Cluster Across Tissues"
        )
    plt.suptitle(title, fontsize=16, weight="bold")
    plt.tight_layout(rect=[0, 0, 0.95, 0.99])
    return fig
