# Introduction

## Overview

This booklet presents a benchmarking analysis comparing two single-nucleus RNA sequencing (snRNA-seq) workflows: **Singulator + FACS** and **Singulator + LeviCell**. The aim is to create a transparent, reproducible, and modular template that others can adapt to benchmark alternative protocols in their own experiments.

This work was carried out as part of the ongoing effort by the **Single Cell Analysis Innovation Lab (SAIL)** to optimize tissue dissociation and nuclei isolation strategies for robust and scalable single-nucleus transcriptomic profiling.

## Why This Comparison?

While both FACS (Fluorescence-Activated Cell Sorting) and LeviCell (microfluidics-based label-free sorting) can enrich for viable nuclei, they differ in principles and potential biases. 

- **FACS** sorts nuclei based on DAPI or other fluorescent labels and can be more selective, but may introduce mechanical stress.
- **LeviCell**, in contrast, uses acoustic waves to separate cells and nuclei without labels or direct contact, which may preserve more fragile populations.

By directly comparing these two protocols using matched input samples and harmonized preprocessing pipelines, we aim to:

- Quantify differences in **cell/nucleus recovery**, **gene detection**, and **transcriptomic complexity**
- Assess potential biases in **cell type composition**
- Identify trade-offs in **throughput**, **cost**, and **technical artifacts**
- Build a **benchmarking template** for future method comparisons

## Data Collection

Matched liver and colon tissue samples were dissociated using the Singulator platform, followed by either FACS or LeviCell enrichment of nuclei. Library preparation and sequencing were performed identically across conditions.

All data were collected at **SAIL**, and the raw and processed data files are stored on the **Iris** file system under:

`/data1/collab002/sail/isabl/datalake/prod/010/collaborators/SAIL/projects/singulator_debris_removal_and/experiments`

## Structure of This Book

This Jupyter Book is organized into modular chapters, each covering a specific stage of the benchmarking pipeline:

- **Introduction**: Current section, providing background and context  
- **Data Preprocessing**: QC, filtering, and normalization  
- **Cell Type Annotation**: Marker-based annotation and scoring  
- **Technical Metrics**: Comparison of read depth, UMI counts, and gene detection  
- **Biological Metrics**: Differences in cell type composition and transcriptomic variance  
- **Doublet Detection and Ambient RNA**: Comparison of noise profiles across protocols  
- **Final Summary**: Recommendations and blueprint for other benchmarking studies

## Reproducibility and Access

All code used in this analysis is publicly available on [GitHub](https://github.com/YOUR_REPO_HERE) with version-controlled Jupyter notebooks. Random seeds are set for reproducibility, and critical decisions are commented inline for clarity. This repository is designed as a blueprint that others can clone and adapt to benchmark their own workflows.
