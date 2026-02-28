# Trajectory analysis of SARS-Cov-2 infection dynamics

Ravindra et al., 2021 performed single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium to identify target cells, alterations in gene expression, and cell state changes (article: `https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001143`), with the data deposited at: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166766`.

The task here is to reproduce the neighbourhood clustering and cell type identification in the referenced paper above (Figures 3A, 3B, 4A, 4B), and also perform pseudotime analysis to order the differentiation of cells. The steps covered in this project include the following:

1. Installation of packages
2. Loading and preprocessing of data  
3. Quality control (QC)  
4. Normalisation and feature selection  
5. Dimensionality reduction (PCA, UMAP)  
6. Clustering  
7. Cell type annotation and visualisations
8. Trajectory inference and pseudotime analysis
9. Biological interpretation of results

## Installation
All necessary Python packages were installed as follows:
```py
!pip install scanpy
!pip install anndata
!pip3 install igraph
!pip install celltypist
!pip install decoupler
!pip install fa2-modified
!pip install louvain
```

## Loading Data

```py
# Import core single cell datasets
import scanpy as sc
import anndata as ad
import numpy as np
```

The raw data (GSE166766_RAW.tar) was downloaded (from `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE166766`) to the local device and uploaded into a Google Drive which was mounted as follows:
```py
# Mounting Google Drive to assess downloaded datasets
from google.colab import drive
drive.mount('/content/drive')
```

Thereafter, the raw data was untarred and re-organised into separate folders according to each disease condition (mock, 1dpi, 2dpi, and 3dpi). This was necessary to be able to treat each disease condition separately to adequately answer the biological questions afterwards.
```py
# Untar (extract) the raw data

!tar -xvf "/content/drive/MyDrive/Trajectory/GSE166766_RAW.tar" -C "/content/drive/MyDrive/Trajectory/"

# !: Tells Python to run a command as a Linux shell command and not Python code.
# tar: 'Tape ARchive', a Linux utility used to create or extract archive files (.tar, .tar.gz, etc.).
# -xvf: '-x' extracts files from the archive; '-v' (verbose mode) lists files as they are extracted; '-f' (file) specifies the archive file to operate on.
# -C: 'Change directory'. Tells tar to extract the files into specified directory.

import os
# Define the base directory
base_dir = "/content/drive/MyDrive/Trajectory"

# Create condition-specific folders
conditions = {
    "mock": "GSM5082289",
    "1dpi": "GSM5082290",
    "2dpi": "GSM5082291",
    "3dpi": "GSM5082292"
}

for cond in conditions:
    os.makedirs(os.path.join(base_dir, cond), exist_ok=True)

# Move files into the correct folders

import shutil                                   # Imports shell utilities for file operations

for cond, gsm in conditions.items():            # Loops over conditions.
    cond_dir = os.path.join(base_dir, cond)

   for fname in os.listdir(base_dir):           # Lists all files in the base directory.
        if fname.startswith(gsm):               # Matches files belonging to a specific GSM.
            src = os.path.join(base_dir, fname) # Defines the source file path.

            # Assign correct 10x filenames
            if "matrix.mtx" in fname:
                new_name = "matrix.mtx.gz"
            elif "features.tsv" in fname:
                new_name = "features.tsv.gz"
            elif "barcodes.tsv" in fname:
                new_name = "barcodes.tsv.gz"
            else:
                continue                       # safety

            dst = os.path.join(base_dir, cond, new_name) # Defines the destination path (with correct renaming).
            shutil.move(src, dst)

# Import and read datasets separately
mock_adata = sc.read_10x_mtx(
    os.path.join(base_dir, "mock"),
    var_names="gene_symbols",
    cache=True
)

day1_adata = sc.read_10x_mtx(
    os.path.join(base_dir, "1dpi"),
    var_names="gene_symbols",
    cache=True
)

day2_adata = sc.read_10x_mtx(
    os.path.join(base_dir, "2dpi"),
    var_names="gene_symbols",
    cache=True
)

day3_adata = sc.read_10x_mtx(
    os.path.join(base_dir, "3dpi"),
    var_names="gene_symbols",
    cache=True
)

mock_adata = sc.read_10x_mtx("/content/drive/MyDrive/Trajectory/mock")
day1_adata = sc.read_10x_mtx("/content/drive/MyDrive/Trajectory/1dpi")
day2_adata = sc.read_10x_mtx("/content/drive/MyDrive/Trajectory/2dpi")
day3_adata = sc.read_10x_mtx("/content/drive/MyDrive/Trajectory/3dpi")
```
```py
# Annotate condition metadata
mock_adata.obs["condition"] = "mock"
day1_adata.obs["condition"] = "1dpi"
day2_adata.obs["condition"] = "2dpi"
day3_adata.obs["condition"] = "3dpi"
```
```py
# Check datasets
print(mock_adata)
print(day1_adata)
print(day2_adata)
print(day3_adata)

mock_adata.obs.head()
day1_adata.obs.head()
day2_adata.obs.head()
day3_adata.obs.head()

mock_adata.var.head()
day1_adata.var.head()
day2_adata.var.head()
day3_adata.var.head()
```
```py
# Make unique
mock_adata.var_names_make_unique()
day1_adata.var_names_make_unique()
day2_adata.var_names_make_unique()
day3_adata.var_names_make_unique()
```

## Quality Control (QC)
The datasets were good quality data from the original authors of the article. However, quality control analysis was performed separately for each disease condition, to confirm the quality of the datasets and to adhere to the usual scRNA-seq analysis pipeline.

```py
# Checking for possible contamination from dying cells, ribosomal transcripts or haemoglobin

# As a rule of thumb:
# Cells with a high proportion of mitochondrial reads (say >10–20%) are likely stressed, apoptotic, or poorly captured.
# Ribosomal transcripts are removed because they represent global transcriptional activity, not cell-type-specific biology.
# Instead of true cell populations, high HB signal often represents ambient RNA contamination from lysed red blood cells.

# mock
mock_adata.var['MT'] = mock_adata.var_names.str.startswith("MT-")
mock_adata.var['RIBO'] = mock_adata.var_names.str.startswith("RPS", "RPL")
mock_adata.var['HB'] = mock_adata.var_names.str.startswith("^HB[^(P)]")

# 1dpi
day1_adata.var['MT'] = day1_adata.var_names.str.startswith("MT-")
day1_adata.var['RIBO'] = day1_adata.var_names.str.startswith("RPS", "RPL")
day1_adata.var['HB'] = day1_adata.var_names.str.startswith("^HB[^(P)]")

# 2dpi
day2_adata.var['MT'] = day2_adata.var_names.str.startswith("MT-")
day2_adata.var['RIBO'] = day2_adata.var_names.str.startswith("RPS", "RPL")
day2_adata.var['HB'] = day2_adata.var_names.str.startswith("^HB[^(P)]")

# 3dpi
day3_adata.var['MT'] = day3_adata.var_names.str.startswith("MT-")
day3_adata.var['RIBO'] = day3_adata.var_names.str.startswith("RPS", "RPL")
day3_adata.var['HB'] = day3_adata.var_names.str.startswith("^HB[^(P)]")

# Calculating the QC metrics

# mock
sc.pp.calculate_qc_metrics(
    mock_adata, qc_vars=["MT", 'RIBO', 'HB'], inplace=True, log1p=True
)

# 1dpi
sc.pp.calculate_qc_metrics(
    day1_adata, qc_vars=["MT", 'RIBO', 'HB'], inplace=True, log1p=True
)

# 2dpi
sc.pp.calculate_qc_metrics(
    day2_adata, qc_vars=["MT", 'RIBO', 'HB'], inplace=True, log1p=True
)

# 3dpi
sc.pp.calculate_qc_metrics(
    day3_adata, qc_vars=["MT", 'RIBO', 'HB'], inplace=True, log1p=True
)
```

Standardising matplotlib for all plots
```py
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (5,4)  # Adjust figure size
plt.rcParams["axes.grid"] = True  # Add grid to plots
plt.rcParams["axes.edgecolor"] = "black" # Set plot border color
plt.rcParams["axes.linewidth"] = 1.5 # Set plot border width
plt.rcParams["axes.facecolor"] = "white" # Set background color
plt.rcParams["axes.labelcolor"] = "black" # Set label color
plt.rcParams["xtick.color"] = "black" # Set x-axis tick color
plt.rcParams["ytick.color"] = "black" # Set y-axis tick color
plt.rcParams["text.color"] = "black" # Set text color
```

Checking for mitochondrial, ribosomal, and haemoglobin genes

```py
# Mock
mt_genes = mock_adata.var[mock_adata.var['MT']]
mt_genes

ribo_genes = mock_adata.var[mock_adata.var['RIBO']]
ribo_genes

hb_genes = mock_adata.var[mock_adata.var['HB']]
hb_genes

# 1dpi
mt_genes = day1_adata.var[day1_adata.var['MT']]
mt_genes

ribo_genes = day1_adata.var[day1_adata.var['RIBO']]
ribo_genes

hb_genes = day1_adata.var[day1_adata.var['HB']]
hb_genes

# 2dpi
mt_genes = day2_adata.var[day2_adata.var['MT']]
mt_genes

ribo_genes = day2_adata.var[day2_adata.var['RIBO']]
ribo_genes

hb_genes = day2_adata.var[day2_adata.var['HB']]
hb_genes

# 3dpi
mt_genes = day3_adata.var[day3_adata.var['MT']]
mt_genes

ribo_genes = day3_adata.var[day3_adata.var['RIBO']]
ribo_genes

hb_genes = day3_adata.var[day3_adata.var['HB']]
hb_genes
```

There are no mitochondrial, ribosomal, and haemoglobin genes. Verify with plots as follows:

```py
# Mock
sc.pl.violin(
    mock_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_MT'],
    jitter=0.4,
    multi_panel=False,
)

sc.pl.violin(
    mock_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_RIBO'],
    jitter=0.4,
    multi_panel=False,
)

sc.pl.violin(
    mock_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_HB'],
    jitter=0.4,
    multi_panel=False,
)

# 1dpi
sc.pl.violin(
    day1_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_MT'],
    jitter=0.4,
    multi_panel=False,
)

sc.pl.violin(
    day1_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_RIBO'],
    jitter=0.4,
    multi_panel=False,
)

sc.pl.violin(
    day1_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_HB'],
    jitter=0.4,
    multi_panel=False,
)

# 2dpi
sc.pl.violin(
    day2_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_MT'],
    jitter=0.4,
    multi_panel=False,
)

sc.pl.violin(
    day2_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_RIBO'],
    jitter=0.4,
    multi_panel=False,
)

sc.pl.violin(
    day2_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_HB'],
    jitter=0.4,
    multi_panel=False,
)

# 3dpi
sc.pl.violin(
    day3_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_MT'],
    jitter=0.4,
    multi_panel=False,
)

sc.pl.violin(
    day3_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_RIBO'],
    jitter=0.4,
    multi_panel=False,
)

sc.pl.violin(
    day3_adata,
    ["n_genes_by_counts", 'total_counts', 'pct_counts_HB'],
    jitter=0.4,
    multi_panel=False,
)




# Visualise each of them in scatter plots.
sc.pl.scatter(mock_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_MT')
sc.pl.scatter(day1_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_MT')
sc.pl.scatter(day2_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_MT')
sc.pl.scatter(day3_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_MT')

sc.pl.scatter(mock_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_RIBO')
sc.pl.scatter(day1_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_RIBO')
sc.pl.scatter(day2_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_RIBO')
sc.pl.scatter(day3_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_RIBO')

sc.pl.scatter(mock_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_HB')
sc.pl.scatter(day1_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_HB')
sc.pl.scatter(day2_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_HB')
sc.pl.scatter(day3_adata, 'total_counts', 'n_genes_by_counts', color='pct_counts_HB')
```
All the steps above confirmed there no MT, RIBO, and HB genes in the dataset.

### Filtering
No necessary for this dataset, since there were no obvious mitochondrial, ribosomal, and haemoglobin genes.

### Doublet detection
```py
sc.pp.scrublet(mock_adata)
sc.pp.scrublet(day1_adata)
sc.pp.scrublet(day2_adata)
sc.pp.scrublet(day3_adata)
```

## Normalisation
Since normalisation adjusts for sequencing depth differences between cells, counts in the dataset were scaled so that each cell can the same total expression level.
```py
# First, save a copy of the data as follows:
bone_marrow_adata.layers["counts"] = bone_marrow_adata.X.copy()

# Normalising to median total counts
sc.pp.normalize_total(bone_marrow_adata)

# Logarithmize the data
sc.pp.log1p(bone_marrow_adata)
```

## Feature selection
```py
# Selecting the top 1000 most variable genes
sc.pp.highly_variable_genes(bone_marrow_adata, n_top_genes=1000)

sc.pl.highly_variable_genes(bone_marrow_adata )        # Left is normalized; right is not.
```

## Dimensionality Reduction (PCA)
**Principal Component Analysis (PCA)** was used to reduce the data complexity and highlight key variation patterns for faster and more robust downstream steps like clustering and visualisation.

```py
sc.tl.pca(bone_marrow_adata)
sc.pl.pca_variance_ratio(bone_marrow_adata, n_pcs=10, log=False)

sc.pl.pca(
    bone_marrow_adata,
    color=["n_genes_by_counts"]
)

sc.pl.pca(
    bone_marrow_adata,
    color=["total_counts"]
)
```

### Nearest Neighbour
The PCA representation of the data matrix was used to compute the neighborhood graph of cells, to cluster the PCA components.
```py
sc.pp.neighbors(bone_marrow_adata)
sc.tl.umap(bone_marrow_adata)

sc.pl.umap(
    bone_marrow_adata,
    color=["n_genes_by_counts"],
    size=8,
)

sc.pl.umap(
    bone_marrow_adata,
    color=["total_counts"],
    size=8,
)
```

## Clustering by communities.
Grouping of cells that show similar expression profiles, for cell type detection.
```py
# Using the igraph implementation and a fixed number of iterations.
sc.tl.leiden(bone_marrow_adata, flavor="igraph", n_iterations=2)

sc.pl.umap(
    bone_marrow_adata,
    color=["leiden"],
    size=8,
)

# Playing around sizing and spacing
sc.pl.umap(
    bone_marrow_adata,
    color=["leiden"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
    ncols = 1
)

# Checking for doublets
sc.pl.umap(
    bone_marrow_adata,
    color=[ "predicted_doublet"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
    ncols = 1
)

sc.pl.umap(
    bone_marrow_adata,
    color=[ "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
    ncols = 1
)

# Further reclustering
sc.tl.leiden(bone_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res0_02", resolution=0.02)
sc.tl.leiden(bone_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(bone_marrow_adata, flavor="igraph", n_iterations=2, key_added="leiden_res2", resolution=2)

sc.pl.umap(
    bone_marrow_adata,
    color=["leiden_res0_02"],
    # increase horizontal space between panels
    wspace=0.5,
    size=15,
    ncols = 1
)

sc.pl.umap(
    bone_marrow_adata,
    color=["leiden_res0_5"],
    # increase horizontal space between panels
    wspace=0.5,
    size=15,
    ncols = 1,
    legend_loc="on data"
)

sc.pl.umap(
    bone_marrow_adata,
    color=["leiden_res2"],
    # increase horizontal space between panels
    wspace=0.5,
    size=15,
    ncols = 1,
    legend_loc="on data"
)
```
Leiden resolution of `0.5` produced relatively sizeable clusters. So, this resolution was used in the downstream analyses.

## Cell Type Annotation
Assigning biological meaning (e.g., cell type or functional state) to each cluster found after Leiden clustering, using **Decoupler**.

**Data correction**:

The dataset was obtained from CZI (Chan Zuckerberg Institute), and uses ensemble gene ids, contrary to what decoupler expects.

To fix this includes running the following codes before importing decoupler:
```py
!wget wget -O result.txt 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Attribute name = "ensembl_gene_id" /><Attribute name = "external_gene_name" /></Dataset></Query>'

# This downloads the table of genes directly from ensemble
```

```py
import pandas as pd

ensembl_var = pd.read_csv('/content/result.txt', header = None)
ensembl_var.columns = ['ensembl_gene_id', 'gene_name']
ensembl_var.head(3)
```
```py
import decoupler as dc

# Query Omnipath and get PanglaoDB
markers = dc.op.resource(name="PanglaoDB", organism="human")

# Keep canonical cell type markers alone
#markers = markers[markers["canonical_marker"]]

# Remove duplicated entries
markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]

# Format because dc only accepts cell_type and genesymbol
markers = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})
markers = markers[["source", "target"]]

markers.head()
```
```py
# Correct target to ensemble
markers = markers.merge(ensembl_var, left_on="target", right_on="gene_name", how="left")
markers = markers.drop(columns=["target"])

# Remove duplicated entries
markers = markers[~markers.duplicated(["source", "ensembl_gene_id"])]

# Format because dc only accepts cell_type and genesymbol
markers = markers.rename(columns={"source": "source", "ensembl_gene_id": "target"})

markers = markers[["source", "target"]]
markers = markers.dropna()

markers.head()
```
```py
# Load the gene expression matrix into dc
dc.mt.ulm(data=bone_marrow_adata,
          net=markers,
          tmin = 3)

# Retrieve the score for each cell type
score = dc.pp.get_obsm(bone_marrow_adata, key="score_ulm")
score

# Preview the data
bone_marrow_adata.obsm["score_ulm"].head()
bone_marrow_adata.obsm["score_ulm"].columns

sc.pl.umap(score, color=["B cells memory", "leiden_res0_5"], cmap="RdBu_r")
```
```py
import seaborn as sns

sc.pl.violin(score, keys=["B cells memory"], groupby="leiden_res0_5", rotation=90)
```
Identifying what each of the 9 clusters mean:
```py
# Rank genes
bone_marrow_adata_rank = dc.tl.rankby_group(score, groupby="leiden_res0_5", reference="rest", method="t-test_overestim_var")
bone_marrow_adata_rank = bone_marrow_adata_rank[bone_marrow_adata_rank["stat"] > 0]
bone_marrow_adata_rank.head()

cluster_annotations = bone_marrow_adata_rank[bone_marrow_adata_rank["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()

cluster_annotations

bone_marrow_adata.obs['cell_type'] = bone_marrow_adata.obs['leiden_res0_5'].map(cluster_annotations)

# Subsetting for multiple genes in the 'source' column
available_genes = set(bone_marrow_adata.var_names)

b_cell_markers = markers[markers['source'].isin(['B cells memory'])]['target']
b_cell_markers = b_cell_markers[b_cell_markers.isin(available_genes)]

nk_cell_markers = markers[markers['source'].isin(['Natural killer T cells'])]['target']
nk_cell_markers = nk_cell_markers[nk_cell_markers.isin(available_genes)]

t_cells_markers = markers[markers['source'].isin(['T cells'])]['target']
t_cells_markers = t_cells_markers[t_cells_markers.isin(available_genes)]

display(b_cell_markers)
display(nk_cell_markers)
display(t_cells_markers)
```

## Visualising cell types in other ways
```py
marker_genes_dict = {
    "B cells": b_cell_markers.head().tolist(),
    "NK cells": nk_cell_markers.head().tolist(),
    "T cells": t_cells_markers.head().tolist()
}

# Dotplot
sc.pl.dotplot(bone_marrow_adata, marker_genes_dict, "cell_type", dendrogram=True)

# Stacked violin plot
sc.pl.stacked_violin(
    bone_marrow_adata, marker_genes_dict, groupby="leiden_res0_5",  dendrogram=True
)

# Matrix plot
sc.pl.matrixplot(
    bone_marrow_adata,
    marker_genes_dict,
    "leiden_res0_5",
    dendrogram=True,
    cmap="Blues",
)

# Heatmap
sc.pl.heatmap(
    bone_marrow_adata, marker_genes_dict, groupby="leiden_res0_5", cmap="viridis", dendrogram=True
)

# Using genome tracks
sc.pl.tracksplot(bone_marrow_adata, marker_genes_dict, groupby="leiden_res0_5", dendrogram=False)
```



## Biological Interpretation: Questions and Answers
### Q1. What cell types did you identify at the different stages of infection?

### Q2. Can you explain why these cell types correlate with COVID-19 infection?

### Q3. Is ACE2 a good marker for tracking COVID-19 infection rate (based on this dataset)?

### Q4. What is the difference between ENO2 and ACE2 as biomarkers in the two studies?

### Q5. Which cell cluster has the highest abundance of ACE2 expression after 3 dpi and what does that mean biologically (INTERPRET VISUALLY)?

## REFERENCES
Ravindra, N.G., Alfajaro, m.M., Gasque, V., Huston, N.C., Wan, H., Szigeti-Buck, K., Yasumoto, Y., Greaney, A.M., Habet, V., Chow, R.D., Chen, J.S., Wei, J., Filler, R.B., Wang, B., Wang, G., Niklason, L.E., Montgomery, R.R., Eisenbarth, S.C., Chen, S., Williams, A., Iwasaki, A., Horvath, T.L., Foxman, E.F., Pierce, R.W., Pyle, A.M., van Dijk, D., Wilen, C.B. (2021). Single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium identifies target cells, alterations in gene expression, and cell state changes. _PLoS Biol_ 19(3): e3001143. https://doi.org/10.1371/journal.pbio.3001143.

