from typing import Dict, Optional, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


# Convert string colormap to dict for categorical legend if needed
def _make_categorical_palette(vals, cmap_name):
    cmap = mpl.colormaps.get_cmap(cmap_name).resampled(len(vals))
    return {val: mpl.colors.to_hex(cmap(i)) for i, val in enumerate(sorted(vals))}


def _resolve_categorical_palette(values, palette):
    if isinstance(palette, str):
        unique_vals = sorted(values.dropna().unique())
        return _make_categorical_palette(unique_vals, palette)
    return palette


def _process_continuous_values(values, clip_low, clip_high, log_scale):
    values = values.clip(clip_low, clip_high)
    return np.log1p(values) if log_scale else values


def _add_colorbar(scatter, ax, vertical, log_scale, clip_low, clip_high, feature):
    location = "bottom" if vertical else "right"
    orientation = "horizontal" if vertical else "vertical"
    pad = 0.15 if vertical else 0.14

    cb = plt.colorbar(
        scatter,
        ax=ax,
        location=location,
        orientation=orientation,
        fraction=0.046,
        pad=pad,
    )

    if log_scale:
        raw_ticks = np.logspace(np.log10(max(1, clip_low)), np.log10(clip_high), num=6)
        cb.set_ticks(np.log1p(raw_ticks))
        cb.set_ticklabels(
            [f"{int(t/1000)}k" if t >= 1000 else str(int(t)) for t in raw_ticks]
        )
    else:
        raw_ticks = np.linspace(clip_low, clip_high, 6)
        cb.set_ticks(raw_ticks)
        cb.set_ticklabels([f"{int(t):,}" for t in raw_ticks])

    cb.ax.set_ylabel(
        feature.replace("log1p_", "").replace("_", " ").title(),
        rotation=0 if vertical else 270,
        labelpad=-160 if vertical else 15,
        loc="top" if vertical else "center",
    )
    return cb


def _compute_clip_bounds(
    values: pd.Series, clip_values: Optional[Tuple[float, float]], log_scale: bool
) -> Tuple[float, float]:
    if clip_values:
        clip_low = values.quantile(clip_values[0])
        if log_scale:
            clip_low = max(1, clip_low)
        clip_high = values.quantile(clip_values[1])
    else:
        clip_low = values.min()
        clip_high = values.max()
    return clip_low, clip_high


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
    legend_bbox_to_anchor: tuple = None,
    fig_height: float = None,
    fig_width: float = None,
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
            len(tissues),
            1,
            figsize=(fig_width or 6, fig_height or 5 * len(tissues)),
            squeeze=False,
        )
    else:
        fig, axes = plt.subplots(
            1,
            len(tissues),
            figsize=(fig_width or 6 * len(tissues), fig_height or 5),
            squeeze=False,
        )

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
        # Shuffle so that the order is not dependent on the input order
        elif shuffle:
            np.random.seed(0)
            shuffled = np.random.permutation(adata.n_obs)
            coords = coords[shuffled]
            values = values.iloc[shuffled]

        is_categorical = (
            values.dtype.name == "category"
            or values.dtype == object
            or isinstance(color_palette, dict)
        )

        if is_categorical:
            # Handle categorical values
            color_palette = _resolve_categorical_palette(values, color_palette)
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
            # Determine global vmin/vmax
            all_values = pd.concat(
                [adata.obs[feature] for adata in combined_by_tissue.values()]
            )
            clip_low, clip_high = _compute_clip_bounds(
                all_values, clip_values, log_scale
            )

            vmin, vmax = (
                (np.log1p(clip_low), np.log1p(clip_high))
                if log_scale
                else (clip_low, clip_high)
            )

            # Handle continuous values
            values = _process_continuous_values(values, clip_low, clip_high, log_scale)
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
                cb = _add_colorbar(
                    scatter,
                    ax,
                    vertical=vertical,
                    log_scale=log_scale,
                    clip_low=clip_low,
                    clip_high=clip_high,
                    feature=feature,
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

        if legend_bbox_to_anchor is None:
            legend_bbox_to_anchor = (1.05, 1.05) if not vertical else (1.15, 0.95)

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
