# Conclusion

## Overview

This booklet provides a **benchmarking template** for comparing two or more single-cell or single-nucleus protocols. It is structured to be transparent, reproducible, and modular—enabling others to adapt it for benchmarking alternative workflows in their own experimental systems.

## What did we learn from performing these analyses?
What do we know about these protocols? How might we expect these priors to influence our results?

> **FILL IN HERE**  
> List metrics and methods we will use for our analysis
> 
**EXAMPLE (DELETE THIS SECTION BEFORE FINALIZING):**  
-  **Scanpy QC metrics** for library size, gene detection
-  **Apoptosis and Mitochondrial** gene set enrichment to identify low quality cells
-  **Downsampling Reads** for quantifying UMI and gene capture 
-  **Scrublet** for doublet analysis
-  **CellTypist** for cell typing (only for tissues that have associated models)
-  **CellBender** for ambient RNA quantification and correction
-  **InferCNVpy** for tumor cell identification

## Why Use these Metrics to Answer our Overarching Question?
What do we hope to gain from these analyses?

> **FILL IN HERE**  
> Describe the goals for each analysis/metric
>
**EXAMPLE (DELETE THIS SECTION BEFORE FINALIZING):**  
-  **Scanpy QC metrics:** If one method captures cells with higher library size/gene detection, then we will have higher quality cells to draw conclusions from
-  **Apoptosis and Mitochondrial:** If one method is harsher with cells, then it might lead to more cells dying and will obscure the biology
-  **Downsampling Reads:** Will allow us to determine if one method is more efficient (more biologocial signal from less reads)
-  **Scrublet:** Higher doublets might indicate poorer dissociation or more technical issues during sorting and/or loading
-  **CellTypist:** Allows us to quantify the cell types we're able to recover using each method
-  **CellBender:** Will let us know which protocol is more likely to contaminate cells with ambient RNA and compare efficiency in debris removal
-  **InferCNVpy:** Tumor cells are usually more fragile, so seeing which method is better for these cells tha tmay be harder to capture will be useful.

## Data Collection

> **FILL IN YOUR EXPERIMENTAL DESIGN HERE**  
> Describe where the sample came from, who performed the experiments, the protocols being used, and anything else about the data collection design that may be relevant to analysis

> **EXAMPLE (DELETE THIS SECTION):**  
- Matched colon tissue samples were collected from the Karuna Ganesh lab (unsure if they are from the same donor)
- Normal liver tissue was also collected from the Ganesh lab but is separate from the tumor data
- The tissues were processed by Maria Bikou at SAIL
- Tissue was dissociated using the Singulator platform, then sorted using FACS or LeviCell
  
## File Locations

> **FILL IN YOUR FILE SYSTEM INFO HERE**  
> Provide a path here for locating raw and processed files.

**EXAMPLE (DELETE):**  
All files are located on the **Iris** system at:  
`/data1/collab002/sail/isabl/datalake/prod/010/collaborators/SAIL/projects/singulator_debris_removal_and/experiments`


## Reproducibility and Access

All code and notebooks used in this benchmark should be version-controlled and openly accessible.

> **FILL IN YOUR GITHUB OR DATA ACCESS INFO HERE**  
> Include links to your GitHub repository

> **EXAMPLE (DELETE):**  
> All code is available at [https://github.com/sail-mskcc/sail-benchmarking-template#](https://github.com/sail-mskcc/sail-benchmarking-template).  

This resource is designed to serve as a **benchmarking blueprint**—simply replace the example text with your own methods, data, and results.

---
