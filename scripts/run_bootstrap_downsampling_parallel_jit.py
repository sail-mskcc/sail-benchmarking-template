import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pprint import pprint
from typing import Dict, List, Tuple

import dill
import numpy as np
import polars as pl

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(CURRENT_DIR, "../data")
ANALYSIS_DIR = os.path.join(DATA_DIR, "analysis")

OUTPUT_DIR = os.path.join(DATA_DIR, "downsampling_bootstraps")
# Sample metadata - in data/metadata.tsv
SAMPLES = {
    "SF_N": ("MB-4027_SF_N", "Normal Colon", "Singulator+FACS"),
    "SL_N": ("MB-4027_SL_N", "Normal Colon", "Singulator+LeviCell"),
    "SF_T": ("MB-4027_SF_T", "Tumor Colon", "Singulator+FACS"),
    "SL_T": ("MB-4027_SL_T", "Tumor Colon", "Singulator+LeviCell"),
    "SF_LN": ("MB-4027_SF_LN", "Normal Liver", "Singulator+FACS"),
    "SL_LN": ("MB-4027_SL_LN", "Normal Liver", "Singulator+LeviCell"),
}


def run_downsampling(
    molecule_infos: Dict[str, pl.DataFrame], target_total_reads: Dict[str, int]
) -> Dict[str, pl.DataFrame]:
    """
    Runs downsampling on a dictionary of molecule_info DataFrames to a specified total number of reads.

    Args:
        molecule_infos (Dict[str, pl.DataFrame]): Dictionary of Polars DataFrames keyed by sample ID.
        target_total_reads (Dict[str, int]): Dictionary of target total reads keyed by sample ID.

    Returns:
        Dict[str, pl.DataFrame]: Downsampled molecule_info DataFrames.
    """
    downsampled_molecule_infos = {}

    for key, df in molecule_infos.items():
        downsampled_df = perform_downsampling(df, target_total_reads[key])
        downsampled_molecule_infos[key] = downsampled_df

    return downsampled_molecule_infos


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

    gene_dict = {sample_id: [] for sample_id in filtered_molecule_infos}

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

            gene_dict[sample_id].append(list(set(downsampled_df["gene_id"].to_list())))

    return results, gene_dict


def _single_bootstrap_to_disk_per_sample(
    sample_key: str,
    file_path: str,
    target_total_reads: int,
    n_steps: int,
    min_reads: int,
    output_dir: str,
    bootstrap_idx: int,
):
    # Load only this sample’s data
    mol_info = pl.read_parquet(file_path)

    # Run stats
    results, downsampled = compute_gene_discovery_downsampled_stats(
        {sample_key: mol_info}, {sample_key: target_total_reads}, n_steps, min_reads
    )

    # Save per-bootstrap per-sample data
    with open(
        os.path.join(
            output_dir,
            sample_key,
            f"{sample_key}_bootstrap_{bootstrap_idx}_downsampled_molecule_info.dill",
        ),
        "wb",
    ) as f:
        dill.dump(downsampled, f)

    with open(
        os.path.join(
            output_dir,
            sample_key,
            f"{sample_key}_bootstrap_{bootstrap_idx}_downsampled_stats.dill",
        ),
        "wb",
    ) as f:
        dill.dump(results, f)

    return results


def run_bootstrap_downsampling_samplewise_parallel(
    sample_file_paths: Dict[str, str],
    target_total_reads: Dict[str, int],
    output_dir: str,
    n_bootstraps=10,
    n_steps=10,
    min_reads=10_000_000,
):
    os.makedirs(output_dir, exist_ok=True)

    futures = []
    with ProcessPoolExecutor(mp_context=get_context("spawn")) as executor:
        for bootstrap_idx in range(n_bootstraps):
            for sample_key, file_path in sample_file_paths.items():
                futures.append(
                    executor.submit(
                        _single_bootstrap_to_disk_per_sample,
                        sample_key,
                        file_path,
                        target_total_reads[sample_key],
                        n_steps,
                        min_reads,
                        output_dir,
                        bootstrap_idx,
                    )
                )

        results = [f.result() for f in futures]

    # Aggregate only summary stats (skip loading downsampled DFs for now)
    summary = defaultdict(lambda: defaultdict(list))
    for res in results:
        for sample_key, stat in res.items():
            for metric, values in stat.items():
                summary[sample_key][metric].append(values)

    summary_stats = defaultdict(dict)
    for sample_key, metrics in summary.items():
        for metric, values in metrics.items():
            arr = np.array(values)
            summary_stats[sample_key][f"{metric}_mean"] = arr.mean(axis=0)
            summary_stats[sample_key][f"{metric}_std"] = arr.std(axis=0)

    return summary_stats


def main():
    # Read in filtered molecule_info DataFrames
    sample_file_paths = {
        key: os.path.join(ANALYSIS_DIR, key, f"{key}_filtered_molecule_info.parquet")
        for key in SAMPLES
    }

    # Define target total reads for each sample
    target_total_reads = {
        "SF_N": 68_000_000,
        "SL_N": 68_000_000,
        "SF_T": 53_000_000,
        "SL_T": 53_000_000,
        "SF_LN": 115_000_000,
        "SL_LN": 115_000_000,
    }

    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Create output subdirectories
    for sample_key in sample_file_paths.keys():
        sample_output_dir = os.path.join(OUTPUT_DIR, sample_key)
        os.makedirs(sample_output_dir, exist_ok=True)

    # Run bootstrap downsampling
    n_bootstraps = 10  # Set to 1 for testing, change to 10 for full run
    n_steps = 10
    min_reads = 10_000_000

    summary_stats = run_bootstrap_downsampling_samplewise_parallel(
        sample_file_paths,
        target_total_reads,
        output_dir=OUTPUT_DIR,
        n_bootstraps=n_bootstraps,
        n_steps=n_steps,
        min_reads=min_reads,
    )

    with open(
        os.path.join(OUTPUT_DIR, "summary_stats_across_bootstraps.dill"), "wb"
    ) as f:
        dill.dump(summary_stats, f)


if __name__ == "__main__":
    main()
