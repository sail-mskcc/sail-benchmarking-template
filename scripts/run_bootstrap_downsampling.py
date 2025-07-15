import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pprint import pprint
from typing import Dict, List

import dill
import numpy as np
import polars as pl

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
# File and directory paths
ANALYSIS_DIR = os.path.join(CURRENT_DIR, "../data/analysis")
# Sample metadata - in data/metadata.tsv
SAMPLES = {
    "SF_N": ("MB-4027_SF_N", "Normal Colon", "Singulator+FACS"),
    "SL_N": ("MB-4027_SL_N", "Normal Colon", "Singulator+LeviCell"),
    "SF_T": ("MB-4027_SF_T", "Tumor Colon", "Singulator+FACS"),
    "SL_T": ("MB-4027_SL_T", "Tumor Colon", "Singulator+LeviCell"),
    "SF_LN": ("MB-4027_SF_LN", "Normal Liver", "Singulator+FACS"),
    "SL_LN": ("MB-4027_SL_LN", "Normal Liver", "Singulator+LeviCell"),
}


def sample_reads(x, n):
    """
    Subsample counts.

    Example:
    # data looks like
    x = [1, 2, 3, 1]
    total_counts = 7

    # draw 4 counts
    sample_indices = [0, 1, 0, 1, 1, 0, 1]  # length 7, keep 1s, drop 0s

    # get indices
    count_indices = [0, 1, 1, 2, 2, 2, 3]  # length 7, used for grouping

    # group by indices and sum
    sample_counts = [0, 1, 2, 1]
    """
    assert n < x.sum(), f"n ({n}) must be smaller than total counts ({x.sum()})"
    total_count = x.sum()
    sample_indices = np.random.binomial(1, n / total_count, total_count)
    count_indices = np.repeat(np.arange(len(x)), x)
    sample_counts = np.bincount(count_indices, weights=sample_indices).astype("int32")
    return sample_counts


def perform_downsampling(
    mol_info_df: pl.DataFrame, target_total_reads: int
) -> pl.DataFrame:
    """
    Downsamples a molecule_info DataFrame to a specified total number of reads.

    Args:
        mol_info_df (pl.DataFrame): DataFrame with a 'reads' column.
        target_total_reads (int): Desired total number of reads after downsampling.

    Returns:
        pl.DataFrame: Filtered and downsampled molecule_info DataFrame.
    """
    read_counts = mol_info_df["reads"].to_numpy()
    downsampled_reads = sample_reads(read_counts, target_total_reads)

    # Keep only molecules with >0 downsampled reads
    nonzero_mask = downsampled_reads > 0
    filtered_df = mol_info_df.filter(pl.Series(nonzero_mask))
    filtered_df = filtered_df.with_columns(
        pl.Series("reads", downsampled_reads[nonzero_mask])
    )

    return filtered_df


def read_in_filtered_molecule_infos():
    """
    Reads in filtered molecule_info DataFrames from parquet files.

    Returns:
        Dict[str, pl.DataFrame]: Dictionary of filtered molecule_info DataFrames keyed by sample ID.
    """
    filtered_molecule_infos = {}
    for key in SAMPLES.keys():
        file_path = os.path.join(
            ANALYSIS_DIR, key, f"{key}_filtered_molecule_info.parquet"
        )
        df = pl.read_parquet(file_path)
        filtered_molecule_infos[key] = df
    return filtered_molecule_infos


def compute_gene_discovery_downsampled_stats(
    filtered_molecule_infos: Dict[str, pl.DataFrame],
    target_total_reads: Dict[str, int],
    n_steps: int = 10,
    min_reads: int = 10_000_000,
) -> Dict[str, Dict[str, list]]:
    """
    Compute gene discovery statistics (median UMIs, genes per cell, and mean reads per cell) by downsampling.

    Args:
        filtered_molecule_infos (Dict[str, pl.DataFrame]): Dict of molecule_info DataFrames keyed by sample ID.
        target_total_reads (Dict[str, int]): Dict of target total reads per sample.
        n_steps (int): Number of evenly spaced read levels to compute between min_reads and max.
        min_reads (int): Minimum number of reads to start downsampling from.

    Returns:
        Dict[str, Dict[str, list]]: Results keyed by sample with read depths and median stats.
    """
    results = {
        sample_id: {
            "total_reads": [],
            "mean_reads": [],
            "median_umis": [],
            "median_genes": [],
        }
        for sample_id in filtered_molecule_infos
    }

    downsampled_df_dict = {sample_id: [] for sample_id in filtered_molecule_infos}

    for sample_id, mol_info in filtered_molecule_infos.items():
        max_reads = target_total_reads[sample_id]
        steps = np.linspace(min_reads, max_reads, n_steps, dtype=int)

        for target in steps:
            # Downsample reads
            downsampled_df = perform_downsampling(mol_info, target)

            # Group by cell once, computing total reads and mean reads
            cell_reads = downsampled_df.group_by("cell").agg(
                pl.sum("reads").alias("reads")
            )

            total_reads = cell_reads["reads"].sum()
            mean_reads = cell_reads["reads"].mean()

            # Compute unique UMIs per (cell, gene)
            umi_counts = downsampled_df.group_by(["cell", "gene_id"]).agg(
                pl.col("umi").n_unique().alias("umi_count")
            )

            # Group again by cell for final summaries
            grouped = umi_counts.group_by("cell").agg(
                [
                    pl.sum("umi_count").alias("umis"),
                    pl.count("gene_id").alias("genes"),
                ]
            )

            median_umis = grouped["umis"].median()
            median_genes = grouped["genes"].median()

            # Append to results
            results[sample_id]["total_reads"].append(float(total_reads))
            results[sample_id]["mean_reads"].append(float(mean_reads))
            results[sample_id]["median_umis"].append(float(median_umis))
            results[sample_id]["median_genes"].append(float(median_genes))

            downsampled_df_dict[sample_id].append(downsampled_df)

    return results, downsampled_df_dict


def run_bootstrap_downsampling(
    filtered_molecule_infos: Dict[str, pl.DataFrame],
    target_total_reads: Dict[str, int],
    n_bootstraps: int = 10,
    n_steps: int = 10,
    min_reads: int = 10_000_000,
) -> Dict[str, Dict[str, List]]:
    """
    Runs bootstrap downsampling on filtered molecule_info DataFrames.

    Args:
        n_bootstraps (int): Number of bootstrap iterations.
        n_steps (int): Number of evenly spaced read levels to compute.
        min_reads (int): Minimum number of reads to start downsampling from.

    Returns:
        Dict[str, Dict[str, List]]: Bootstrap results and downsampled DataFrames.
    """

    # Run 10 bootstraps and store each replicate
    bootstrap_results = defaultdict(lambda: defaultdict(list))
    downsampled_df_bootstrap_list = []

    for i in range(n_bootstraps):
        bootstrap_result, downsampled_df_dict = (
            compute_gene_discovery_downsampled_stats(
                filtered_molecule_infos=filtered_molecule_infos,
                target_total_reads=target_total_reads,
                n_steps=n_steps,
                min_reads=min_reads,
            )
        )
        downsampled_df_bootstrap_list.append(downsampled_df_dict)
        for sample_id in bootstrap_result:
            for key in ["total_reads", "mean_reads", "median_umis", "median_genes"]:
                bootstrap_results[sample_id][key].append(
                    bootstrap_result[sample_id][key]
                )

    # Summarize into plot-ready format
    summary_stats = defaultdict(dict)
    for sample_id, stats in bootstrap_results.items():
        # Each key is a list of n_bootstraps × n_steps
        total_reads_array = np.array(
            stats["total_reads"]
        )  # shape: (n_bootstraps, n_steps)
        mean_reads_array = np.array(
            stats["mean_reads"]
        )  # shape: (n_bootstraps, n_steps)
        median_umis_array = np.array(stats["median_umis"])
        median_genes_array = np.array(stats["median_genes"])

        # Mean and std across bootstraps
        summary_stats[sample_id]["total_reads_mean"] = total_reads_array.mean(axis=0)
        summary_stats[sample_id]["total_reads_std"] = total_reads_array.std(axis=0)

        summary_stats[sample_id]["mean_reads_mean"] = mean_reads_array.mean(axis=0)
        summary_stats[sample_id]["mean_reads_std"] = mean_reads_array.std(axis=0)

        summary_stats[sample_id]["median_umis_mean"] = median_umis_array.mean(axis=0)
        summary_stats[sample_id]["median_umis_std"] = median_umis_array.std(axis=0)

        summary_stats[sample_id]["median_genes_mean"] = median_genes_array.mean(axis=0)
        summary_stats[sample_id]["median_genes_std"] = median_genes_array.std(axis=0)

    return bootstrap_results, downsampled_df_dict, summary_stats


def main():
    # Read in filtered molecule_info DataFrames
    filtered_molecule_infos = read_in_filtered_molecule_infos()

    # Define target total reads for each sample
    target_total_reads = {
        "SF_N": 68_000_000,
        "SL_N": 68_000_000,
        "SF_T": 53_000_000,
        "SL_T": 53_000_000,
        "SF_LN": 115_000_000,
        "SL_LN": 115_000_000,
    }

    # Run bootstrap downsampling
    n_bootstraps = 10  # Set to 1 for testing, change to 10 for full run
    n_steps = 10
    min_reads = 10_000_000

    bootstrap_results, downsampled_df_bootstrap_list, summary_stats = (
        run_bootstrap_downsampling(
            filtered_molecule_infos=filtered_molecule_infos,
            target_total_reads=target_total_reads,
            n_bootstraps=n_bootstraps,
            n_steps=n_steps,
            min_reads=min_reads,
        )
    )

    pprint(f"{bootstrap_results=}")
    pprint(f"{downsampled_df_bootstrap_list=}")
    pprint(f"{summary_stats=}")

    # Structure - list with len(n_bootstraps), each element is a dict with sample_id keys, each sample_id has a dict with keys and a list of n_steps dicts
    with open(
        os.path.join(
            "/data1/collab002/sail/projects/ongoing/benchmarks/sail-benchmarking-template/data",
            "other_script_run",
            "downsampling_molecule_infos.dill",
        ),
        "wb",
    ) as f:
        dill.dump(downsampled_df_bootstrap_list, f)

    with open(
        os.path.join(
            "/data1/collab002/sail/projects/ongoing/benchmarks/sail-benchmarking-template/data",
            "other_script_run",
            "downsampling_bootstrap_results.dill",
        ),
        "wb",
    ) as f:
        dill.dump(bootstrap_results, f)

    with open(
        os.path.join(
            "/data1/collab002/sail/projects/ongoing/benchmarks/sail-benchmarking-template/data",
            "other_script_run",
            "downsampling_summary_stats.dill",
        ),
        "wb",
    ) as f:
        dill.dump(summary_stats, f)


if __name__ == "__main__":
    main()
