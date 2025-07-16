import os
from pathlib import Path

# Current directory - where this file is in
CURRENT_DIR = Path(__file__).resolve().parent

# Data Directories

# Original data directory
ORIGINAL_DATA_ROOT_DIR  = "/data1/collab002/sail/projects/ongoing/benchmarks/benchmark_facs_levicell/data/read_only"

# Local data directories
DATA_DIR = CURRENT_DIR / ".." / "data"
ADATA_DIR = DATA_DIR / "adatas"
GENE_SET_INFO_DIR = DATA_DIR / "gene_set_info"

# Figures Directory
FIGURES_DIR = CURRENT_DIR / "figures"

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

TISSUE_ORDER = ("Normal Colon", "Tumor Colon", "Normal Liver")

# File Names
CELL_RANGER_FILTERED_FEATURE_MATRIX_FILE_NAME = "filtered_feature_bc_matrix.h5"
CELL_RANGER_RAW_FEATURE_MATRIX_FILE_NAME = "raw_feature_bc_matrix.h5"

CELL_RANGER_FILTERED_BARCODES_FILE_NAME = ""