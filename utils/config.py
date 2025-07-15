import os
from pathlib import Path

# Current directory - where this file is in
CURRENT_DIR = Path(__file__).resolve().parent

# Data Directories

# Original data directory
# ORIGINAL_DATA_ROOT_DIR = "/data1/collab002/sail/isabl/datalake/prod/010/collaborators/SAIL/projects/singulator_debris_removal_and/experiments"
ORIGINAL_DATA_ROOT_DIR  = "/data1/collab002/sail/projects/ongoing/benchmarks/benchmark_facs_levicell/data/read_only"

# Local data directories
DATA_DIR = CURRENT_DIR / ".." / "data"
ADATA_DIR = DATA_DIR / "adatas"

# Sample Metadata
SAMPLE_IDENTIFIER = "MB-4027"
SAMPLES_METADATA = {
    "SF_N": {"tissue" : "Normal Colon", "protocol" : "Singulator+FACS"},
    "SL_N": { "tissue" : "Normal Colon", "protocol" : "Singulator+LeviCell"},
    "SF_T": {"tissue" : "Tumor Colon", "protocol" : "Singulator+FACS"},
    "SL_T": {"tissue" : "Tumor Colon", "protocol" : "Singulator+LeviCell"},
    "SF_LN": {"tissue" : "Normal Liver", "protocol" : "Singulator+FACS"},
    "SL_LN": {"tissue" : "Normal Liver", "protocol" : "Singulator+LeviCell"},
}

# Color palette for plotting
PROTOCOL_COLOR_PALETTE = {
    "Singulator+FACS": "#AEC6CF",
    "Singulator+LeviCell": "#FFDAB9",
}

# File Names
CELL_RANGER_FILTERED_FEATURE_MATRIX_FILE_NAME = "filtered_feature_bc_matrix.h5"
CELL_RANGER_RAW_FEATURE_MATRIX_FILE_NAME = "raw_feature_bc_matrix.h5"

# # Original data root directory
# ORIGINAL_DATA_ROOT_DIR = "`/data1/collab002/sail/isabl/datalake/prod/010/collaborators/SAIL/projects/singulator_debris_removal_and/experiments"


# GENE_SETS_DIR = DATA_DIR / "gene_sets"

# CELL_RANGER_FILTERED_FEATURE_MATRIX_FILE_NAME = "filtered_feature_bc_matrix.h5"
# CELL_RANGER_RAW_FEATURE_MATRIX_FILE_NAME = "raw_feature_bc_matrix.h5"


# ANALYSIS_DIR = DATA_DIR / "analysis"
# INPUT_ADATA_DIR = ANALYSIS_DIR / "adatas" / "adatas_X_cellbender_raw_filtered"
# OUTPUT_ADATA_DIR = (
#     ANALYSIS_DIR / "adatas" / "adatas_X_cellbender_raw_filtered_with_gene_metrics"
# )

# # Make sure output adata dir directory exists
# os.makedirs(OUTPUT_ADATA_DIR, exist_ok=True)

# # Sample metadata - in data/metadata.tsv
# samples = {
#     "SF_N": ("MB-4027_SF_N", "Normal Colon", "Singulator+FACS"),
#     "SL_N": ("MB-4027_SL_N", "Normal Colon", "Singulator+LeviCell"),
#     "SF_T": ("MB-4027_SF_T", "Tumor Colon", "Singulator+FACS"),
#     "SL_T": ("MB-4027_SL_T", "Tumor Colon", "Singulator+LeviCell"),
#     "SF_LN": ("MB-4027_SF_LN", "Normal Liver", "Singulator+FACS"),
#     "SL_LN": ("MB-4027_SL_LN", "Normal Liver", "Singulator+LeviCell"),
# }

# # Color palette for plotting
# protocol_color_palette = {
#     "Singulator+FACS": "#AEC6CF",
#     "Singulator+LeviCell": "#FFDAB9",
# }
