#!/bin/bash


CURRENT_DIR=$(realpath "$(dirname "$0")")
ANALYSIS_DIR=$(realpath "$CURRENT_DIR/../data/analysis")

module load cuda/12.0

if ! mamba env list | grep -q '^cellbender_env'; then
    echo "Environment 'cellbender_env' does not exist."
    mamba create -n cellbender_env python=3.7 -y
    mamba activate cellbender_env

    mamba install -c conda-forge pytables -y
    pip install lxml_html_clean
    pip3 install torch torchvision torchaudio
    pip install cellbender

    mamba deactivate
else
    echo "Environment 'cellbender_env' already exists."
fi



for subdir in $ANALYSIS_DIR/*; do
    # Skip if it's not a directory or if it's 'adatas'
    [ -d "$subdir" ] || continue
    [ "$(basename "$subdir")" = "adatas" ] && continue

    # Construct path to cell_bender directoryc
    cell_bender_dir="$subdir/cell_bender"
    script_name="$(basename "$subdir")_cell_bender_gpu.sh"
    script_path="$cell_bender_dir/$script_name"
    chmod +x "$script_path"

    # Check if script exists before submitting
    if [[ -f "$script_path" ]]; then
        cd "$cell_bender_dir"
        sbatch "$script_name"
        cd - > /dev/null
    else
        echo "Script not found: $script_path"
    fi
done

# echo "All jobs submitted. Check the cell_bender directories for job scripts."