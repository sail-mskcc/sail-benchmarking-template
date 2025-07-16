from typing import Callable, Dict, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from .common_plotting_utils import adaptive_formatter


def plot_scalar_metric(
    data_dict: dict,
    metadata: dict,
    metric_func: Callable,
    metric_label: str,
    title: str = None,
    annotate: bool = True,
    protocol_color_palette: dict = None,
    adaptive_formatter: Callable = adaptive_formatter,
    tissue_order: list | Tuple = None,
):
    """
    Generic function to plot scalar metrics computed on each sample's data (AnnData or DataFrame),
    grouped by tissue and protocol.

    Args:
        data_dict: dict of sample data (AnnData or DataFrames), keyed by sample ID
        metadata: dict mapping sample ID to {"tissue": ..., "protocol": ...}
        metric_func: function that computes scalar metric for each sample's data
        metric_label: string for Y-axis label and plot title
        protocol_color_palette: dict mapping protocol to colors
        adaptive_formatter: function for formatting y-axis tick labels
        tissue_order: optional list of tissue names to enforce subplot order
    """
    # Assemble data
    rows = []
    for key, sample in data_dict.items():
        tissue = metadata[key]["tissue"]
        protocol = metadata[key]["protocol"]
        value = metric_func(sample)
        rows.append({"Tissue": tissue, "Protocol": protocol, "Value": value})
    df = pd.DataFrame(rows)

    # Determine tissue order
    if tissue_order is None:
        tissue_order = df["Tissue"].unique().tolist()

    # Filter DataFrame to only include tissues in the specified order (preserving order)
    df["Tissue"] = pd.Categorical(df["Tissue"], categories=tissue_order, ordered=True)
    df = df.sort_values("Tissue")

    unique_protocols = df["Protocol"].unique()

    # Create figure
    fig, axes = plt.subplots(
        1, len(tissue_order), figsize=(5 * len(tissue_order), 5), sharey=True
    )
    if len(tissue_order) == 1:
        axes = [axes]

    # Create color palette if not provided
    if protocol_color_palette is None:
        protocol_color_palette = sns.color_palette("husl", len(unique_protocols))
        protocol_color_palette = dict(zip(unique_protocols, protocol_color_palette))

    # Plot each tissue
    for i, tissue in enumerate(tissue_order):
        ax = axes[i]
        df_tissue = df[df["Tissue"] == tissue]
        sns.barplot(
            data=df_tissue,
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

        if annotate:
            for p in ax.patches:
                val = p.get_height()
                ax.annotate(
                    f"{int(val):,}",
                    (p.get_x() + p.get_width() / 2.0, val),
                    ha="center",
                    va="bottom",
                    fontsize=10,
                    weight="bold",
                )

    # Add title and legend
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



def plot_adata_metric_histogram(
    adata_dict: Dict[str, sc.AnnData],
    sample_metadata: Dict[str, tuple],
    field: str,
    axis: str = "obs",
    bins: int = 100,
    row_label: str = "Cell",
    x_label: str = None,
    log_x: bool = True,
    log_y: bool = False,
    title: str = None,
    proportion: bool = False,
    protocol_color_palette: Dict[str, str] = None,
    nonzero: bool = True,
    tissue_order: list | Tuple = None,
) -> None:
    """
    Plots the distribution of a field (column) from each AnnData object, grouped by tissue and protocol.

    Args:
        adata_dict: Dictionary of AnnData objects keyed by sample ID (e.g., 'sf_ln').
        sample_metadata: Dictionary mapping sample ID to (tissue, protocol) tuples.
        field: Column name to plot from .obs or .var (e.g., 'total_counts').
        axis: Whether to pull from .obs or .var (default: 'obs').
        bins: Number of histogram bins.
        x_label: Custom label for the x-axis.
        log_x: Whether to apply log1p transformation to the x-axis values (default: True).
        log_y: Whether to apply log scale to the y-axis (default: False).
        title: Title of the entire figure.
        row_label: Label for the y-axis.
        protocol_color_palette: Optional dictionary mapping protocol names to colors.
        proportion: If True, plot proportions instead of counts.
        tissue_order: optional list of tissue names to enforce subplot order
    """
    # Gather data
    rows = []
    for key, adata in adata_dict.items():
        tissue = sample_metadata[key]["tissue"]
        protocol = sample_metadata[key]["protocol"]
        values = getattr(adata, axis)[field]
        if log_x:
            values = np.log1p(values)
        for val in values:
            rows.append({"Tissue": tissue, "Protocol": protocol, "Value": val})

    df = pd.DataFrame(rows)

    # Determine tissue order
    if tissue_order is None:
        tissue_order = df["Tissue"].unique().tolist()

    # Filter DataFrame to only include tissues in the specified order (preserving order)
    df["Tissue"] = pd.Categorical(df["Tissue"], categories=tissue_order, ordered=True)
    df = df.sort_values("Tissue")

    # Collect unique tissues and protocols
    unique_tissues = df["Tissue"].unique()
    unique_protocols = df["Protocol"].unique()

    # Create color palette if not provided
    if protocol_color_palette is None:
        protocol_color_palette = sns.color_palette("husl", len(unique_protocols))
        protocol_color_palette = dict(zip(unique_protocols, protocol_color_palette))

    # Plot
    sns.set_theme(style="white", font_scale=1.1)
    fig, axes = plt.subplots(
        1, len(unique_tissues), figsize=(5 * len(unique_tissues), 5), sharey=True
    )
    if len(unique_tissues) == 1:
        axes = [axes]

    # Determine stat type
    stat_type = "probability" if proportion else "count"

    if nonzero:
        df = df[df["Value"] > 0]

    # Plot each tissue
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
            [getattr(adata_dict[key], axis)[field].to_numpy() for key in adata_dict]
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
            [getattr(adata_dict[key], axis)[field].to_numpy() for key in adata_dict]
        )
        vmin, vmax = all_values.min(), all_values.max()
        for ax in axes:
            ax.set_xlim(vmin, vmax)

    #  Add Title
    if not title:
        title = f"{field} by Protocol Across Tissues"
    plt.suptitle(title, fontsize=16, weight="bold")

    # Add legend
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


def plot_adata_metric_violin(
    adata_dict: Dict[str, sc.AnnData],
    sample_metadata: Dict[str, tuple],
    field: str,
    axis: str = "obs",
    log_y: bool = True,
    y_label: str = "Value",
    title: str = None,
    protocol_color_palette: Dict[str, str] = None,
    tissue_order: list | Tuple = None,
) -> None:
    """
    Plots a violin distribution of a given .obs or .var field from each AnnData object,

    Args:
        adata_dict: Dictionary of AnnData objects keyed by sample ID (e.g., 'sf_ln').
        sample_metadata: Dictionary mapping sample ID to (tissue, protocol) tuples.
        field: Column name to plot from .obs or .var (e.g., 'total_counts').
        axis: Whether to pull from .obs or .var (default: 'obs').
        log_y: Whether to apply log1p transformation to the y-axis values (default: True).
        title: Title of the entire figure.
        y_label: Label for the y-axis.
        protocol_color_palette: Optional dictionary mapping protocol names to colors.
        tissue_order: optional list of tissue names to enforce subplot order
    """
    # Assemble data
    rows = []
    for key, adata in adata_dict.items():
        tissue = sample_metadata[key]["tissue"]
        protocol = sample_metadata[key]["protocol"]
        values = getattr(adata, axis)[field]
        if log_y:
            values = np.log1p(values)
        for val in values:
            rows.append({"Tissue": tissue, "Protocol": protocol, "Value": val})
    df = pd.DataFrame(rows)

    # Determine tissue order
    if tissue_order is None:
        tissue_order = df["Tissue"].unique().tolist()

    # Filter DataFrame to only include tissues in the specified order (preserving order)
    df["Tissue"] = pd.Categorical(df["Tissue"], categories=tissue_order, ordered=True)
    df = df.sort_values("Tissue")

    # Get unique tissues and protocols
    unique_tissues = df["Tissue"].unique()
    unique_protocols = df["Protocol"].unique()

    # Create color palette if not provided
    if protocol_color_palette is None:
        protocol_color_palette = sns.color_palette("husl", len(unique_protocols))
        protocol_color_palette = dict(zip(unique_protocols, protocol_color_palette))

    # Setup plot grid
    sns.set_theme(style="white", font_scale=1.1)
    fig, axes = plt.subplots(
        1, len(unique_tissues), figsize=(5 * len(unique_tissues), 5), sharey=True
    )
    if len(unique_tissues) == 1:
        axes = [axes]

    # Plot each tissue
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

    # Add title
    if not title:
        title = f"{field} by Protocol Across Tissues"
    plt.suptitle(title, fontsize=16, weight="bold")

    # Custom y-axis ticks if log-transformed
    if log_y:
        raw_values = np.concatenate(
            [getattr(adata_dict[key], axis)[field].to_numpy() for key in adata_dict]
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


def plot_obs_scatter_by_tissue(
    adata_dict: Dict[str, sc.AnnData],
    sample_metadata: Dict[str, tuple],
    x_field: str,
    y_field: str,
    log_x: bool = True,
    log_y: bool = False,
    title: str = "Scatter Plot of Obs Fields by Tissue",
    protocol_color_palette: Dict[str, str] = None,
    shuffle: bool = True,
    tissue_order: list | Tuple = None,
) -> plt.Figure:
    """
    Plots scatter plot of two obs fields for each AnnData object, grouped by tissue.

    Args:
        adata_dict: Dictionary of AnnData objects.
        sample_metadata: Dictionary mapping keys in adata_dict to dict with 'Tissue' and 'Protocol'.
        x_field: Field in adata.obs for x-axis.
        y_field: Field in adata.obs for y-axis.
        log_x: Apply log1p to x-axis values.
        log_y: Apply log1p to y-axis values.
        title: Plot title.
        protocol_color_palette: Optional dictionary for protocol-specific colors.
        shuffle: If True, shuffle the data points to avoid overlap in scatter plot.
        tissue_order: Optional list of tissue names to enforce subplot order.

    Returns:
        matplotlib Figure object.
    """
    rows = []
    for key, adata in adata_dict.items():
        tissue, protocol = sample_metadata[key]["tissue"], sample_metadata[key]["protocol"]
        x = adata.obs[x_field]
        y = adata.obs[y_field]
        if log_x:
            x = np.log1p(x)
        if log_y:
            y = np.log1p(y)
        for xi, yi in zip(x, y):
            rows.append(
                {"Tissue": tissue, "Protocol": protocol, x_field: xi, y_field: yi}
            )
    df = pd.DataFrame(rows)

    if shuffle:
        df = df.sample(frac=1, random_state=0).reset_index(drop=True)

    # Determine tissue order
    if tissue_order is None:
        tissue_order = df["Tissue"].unique().tolist()

    # Filter DataFrame to only include tissues in the specified order (preserving order)
    df["Tissue"] = pd.Categorical(df["Tissue"], categories=tissue_order, ordered=True)
    df = df.sort_values("Tissue")

    unique_tissues = df["Tissue"].unique()
    sns.set_theme(style="white", font_scale=1.1)
    fig, axes = plt.subplots(
        1, len(unique_tissues), figsize=(5 * len(unique_tissues), 5), sharey=True
    )
    if len(unique_tissues) == 1:
        axes = [axes]

    for i, tissue in enumerate(unique_tissues):
        ax = axes[i]
        df_subset = df[df["Tissue"] == tissue]
        if shuffle:
            df_subset = df_subset.sample(frac=1, random_state=0).reset_index(drop=True)
        sns.scatterplot(
            data=df_subset,
            x=x_field,
            y=y_field,
            hue="Protocol",
            palette=protocol_color_palette,
            alpha=0.8,
            size=0.8,
            ax=ax,
            legend=False,
            linewidth=0,
        )
        ax.set_title(tissue, fontsize=13, weight="bold")
        ax.set_xlabel(f"log1p({x_field})" if log_x else x_field)
        if i == 0:
            ax.set_ylabel(f"log1p({y_field})" if log_y else y_field)
        else:
            ax.set_ylabel("")

    plt.suptitle(title, fontsize=16, weight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig