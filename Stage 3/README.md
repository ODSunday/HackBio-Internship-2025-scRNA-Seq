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
# mock
mock_adata.layers["counts"] = mock_adata.X.copy()   # Saves a copy of the data
sc.pp.normalize_total(mock_adata)                   # Normalises to median total counts
sc.pp.log1p(mock_adata)                             # Logarithmizes the data

# 1dpi
day1_adata.layers["counts"] = day1_adata.X.copy()
sc.pp.normalize_total(day1_adata)
sc.pp.log1p(day1_adata)

# 2dpi
day2_adata.layers["counts"] = day2_adata.X.copy()
sc.pp.normalize_total(day2_adata)
sc.pp.log1p(day2_adata)

# 3dpi
day3_adata.layers["counts"] = day3_adata.X.copy()
sc.pp.normalize_total(day3_adata)
sc.pp.log1p(day3_adata)
```

## Feature selection
```py
# Selecting the top 1000 most variable genes
# mock
sc.pp.highly_variable_genes(mock_adata, n_top_genes=1000)
sc.pl.highly_variable_genes(mock_adata)

# 1dpi
sc.pp.highly_variable_genes(day1_adata, n_top_genes=1000)
sc.pl.highly_variable_genes(day1_adata)

# 2dpi
sc.pp.highly_variable_genes(day2_adata, n_top_genes=1000)
sc.pl.highly_variable_genes(day2_adata)

# 3dpi
sc.pp.highly_variable_genes(day3_adata, n_top_genes=1000)
sc.pl.highly_variable_genes(day3_adata)
```

## Dimensionality Reduction (PCA)
**Principal Component Analysis (PCA)** was used to reduce the data complexity and highlight key variation patterns for faster and more robust downstream steps like clustering and visualisation.

```py
# mock
sc.tl.pca(mock_adata)
sc.pl.pca_variance_ratio(mock_adata, n_pcs=10, log=False)

sc.pl.pca(mock_adata,
          color=["condition"],
          cmap="coolwarm")

# 1dpi
sc.tl.pca(day1_adata)
sc.pl.pca_variance_ratio(day1_adata, n_pcs=10, log=False)

sc.pl.pca(day1_adata,
          color=["condition"],
          cmap="coolwarm")

# 2dpi
sc.tl.pca(day2_adata)
sc.pl.pca_variance_ratio(day2_adata, n_pcs=10, log=False)

sc.pl.pca(day2_adata,
          color=["condition"],
          cmap="coolwarm")

# 3dpi
sc.tl.pca(day3_adata)
sc.pl.pca_variance_ratio(day3_adata, n_pcs=10, log=False)

sc.pl.pca(day3_adata,
          color=["condition"],
          cmap="coolwarm")
```
## UMAP
```py
# mock
sc.pp.neighbors(mock_adata)
sc.tl.umap(mock_adata)
mock_adata                       # Views data info.

# 1dpi
sc.pp.neighbors(day1_adata)
sc.tl.umap(day1_adata)
day1_adata

# 2dpi
sc.pp.neighbors(day2_adata)
sc.tl.umap(day2_adata)
day2_adata

# 3dpi
sc.pp.neighbors(day3_adata)
sc.tl.umap(day3_adata)
day3_adata

# Plotting the UMAP by disease codition, ACE2, and ENO2.
# mock
sc.pl.umap(
    mock_adata,
    color=["condition", "ACE2", "ENO2"],
    size=8,
)

# 1dpi
sc.pl.umap(
    day1_adata,
    color=["condition", "ACE2", "ENO2"],
    size=8,
)

# 2dpi
sc.pl.umap(
    day2_adata,
    color=["condition", "ACE2", "ENO2"],
    size=8,
)

# 3dpi
sc.pl.umap(
    day3_adata,
    color=["condition", "ACE2", "ENO2"],
    size=8,
)
```
ACE2 not a good marker for infection, probably.

ENO2 absent in infection states. Absence of ENO2 probably indicates the prsence of infection.

## Clustering

```py
# mock
sc.tl.leiden(mock_adata, flavor="igraph", n_iterations=2, key_added="leiden_res_", resolution=0.25)

# 1dpi
sc.tl.leiden(day1_adata, flavor="igraph", n_iterations=2, key_added="leiden_res_", resolution=0.25)

# 2dpi
sc.tl.leiden(day2_adata, flavor="igraph", n_iterations=2, key_added="leiden_res_", resolution=0.25)

# 3dpi
sc.tl.leiden(day3_adata, flavor="igraph", n_iterations=2, key_added="leiden_res_", resolution=0.25)

# mock
sc.pl.umap(
    mock_adata,
    color=["leiden_res_"],
    size=16,
    title='mock'
)

# 1dpi
sc.pl.umap(
    day1_adata,
    color=["leiden_res_"],
    size=16,
    title='1dpi'
)

# 2dpi
sc.pl.umap(
    day2_adata,
    color=["leiden_res_"],
    size=16,
    title='2dpi'
)

# 3dpi
sc.pl.umap(
    day3_adata,
    color=["leiden_res_"],
    size=16,
    title='3dpi'
)

# ACE2 and ENO2

# mock
sc.pl.umap(
    mock_adata,
    color=["leiden_res_", "ACE2", "ENO2"],
    size=8,
    title='mock'
)

# 1dpi
sc.pl.umap(
    day1_adata,
    color=["leiden_res_", "ACE2", "ENO2"],
    size=8,
    title='1dpi'
)

# 2dpi
sc.pl.umap(
    day2_adata,
    color=["leiden_res_", "ACE2", "ENO2"],
    size=8,
    title='2dpi'
)

# 3dpi
sc.pl.umap(
    day3_adata,
    color=["leiden_res_", "ACE2", "ENO2"],
    size=8,
    title='3dpi'
)

# TMPRSS4

# mock
sc.pl.umap(
    mock_adata,
    color=["leiden_res_", "TMPRSS4"],
    size=8,
    title='mock'
)

# 1dpi
sc.pl.umap(
    day1_adata,
    color=["leiden_res_", "TMPRSS4"],
    size=8,
    title='1dpi'
)

# 2dpi
sc.pl.umap(
    day2_adata,
    color=["leiden_res_", "TMPRSS4"],
    size=8,
    title='2dpi'
)

# 3dpi
sc.pl.umap(
    day3_adata,
    color=["leiden_res_", "TMPRSS4"],
    size=8,
    title='3dpi'
)
```

## Cell Type Annotation
Assigning biological meaning (e.g., cell type or functional state) to each cluster found after Leiden clustering, using **Decoupler**.

```py
import decoupler as dc

# Subset to specific organ to get relevant cell types.

# Query Omnipath and get PanglaoDB
markers = dc.op.resource(name="PanglaoDB", organism="human")

# Inspect 'markers'
markers.shape
markers.head()

# Select a specific organ
markers['organ'].unique()

# Keep canonical cell type markers alone
markers = markers[markers["organ"] == "Lungs"]
markers.shape

# Remove duplicated entries
markers = markers[~markers.duplicated(["cell_type", "genesymbol"])]

# Format because dc only accepts cell_type and genesymbol
markers = markers.rename(columns={"cell_type": "source", "genesymbol": "target"})
markers = markers[["source", "target"]]

markers.head()

mock_adata.var_names
day1_adata.var_names
day2_adata.var_names
day3_adata.var_names

# Load the gene expression matrix into dc
dc.mt.ulm(data=mock_adata,
          net=markers,
          tmin = 3)

dc.mt.ulm(data=day1_adata,
          net=markers,
          tmin = 3)

dc.mt.ulm(data=day2_adata,
          net=markers,
          tmin = 3)

dc.mt.ulm(data=day3_adata,
          net=markers,
          tmin = 3)

# Retrieve the score for each cell type
score_mock = dc.pp.get_obsm(mock_adata, key="score_ulm")
score_1 = dc.pp.get_obsm(day1_adata, key="score_ulm")
score_2 = dc.pp.get_obsm(day2_adata, key="score_ulm")
score_3 = dc.pp.get_obsm(day3_adata, key="score_ulm")

# Preview the data
mock_adata.obsm["score_ulm"].head(3)
day1_adata.obsm["score_ulm"].head(3)
day2_adata.obsm["score_ulm"].head(3)
day3_adata.obsm["score_ulm"].head(3)

mock_adata.obsm["score_ulm"].columns
day1_adata.obsm["score_ulm"].columns
day2_adata.obsm["score_ulm"].columns
day3_adata.obsm["score_ulm"].columns
```

### Rank genes

```py
# mock
mock_adata_gene_rank = dc.tl.rankby_group(score_mock, groupby="leiden_res_", reference="rest", method="t-test_overestim_var")
mock_adata_gene_rank = mock_adata_gene_rank[mock_adata_gene_rank["stat"] > 0]
mock_adata_gene_rank.head(5)

# 1dpi
day1_adata_gene_rank = dc.tl.rankby_group(score_1, groupby="leiden_res_", reference="rest", method="t-test_overestim_var")
day1_adata_gene_rank = day1_adata_gene_rank[day1_adata_gene_rank["stat"] > 0]
day1_adata_gene_rank.head(5)

# 2dpi
day2_adata_gene_rank = dc.tl.rankby_group(score_2, groupby="leiden_res_", reference="rest", method="t-test_overestim_var")
day2_adata_gene_rank = day2_adata_gene_rank[day2_adata_gene_rank["stat"] > 0]
day2_adata_gene_rank.head(5)

# 3dpi
day3_adata_gene_rank = dc.tl.rankby_group(score_3, groupby="leiden_res_", reference="rest", method="t-test_overestim_var")
day3_adata_gene_rank = day3_adata_gene_rank[day3_adata_gene_rank["stat"] > 0]
day3_adata_gene_rank.head(5)
```

### Identify cell types

```py
# mock
top_cell_type_per_group_mock = mock_adata_gene_rank.groupby('group')['name'].apply(lambda x: x.head(1))
display(top_cell_type_per_group.to_dict())

# 1dpi
top_cell_type_per_group_1 = day1_adata_gene_rank.groupby('group')['name'].apply(lambda x: x.head(1))
display(top_cell_type_per_group.to_dict())

# 2dpi
top_cell_type_per_group_2 = day2_adata_gene_rank.groupby('group')['name'].apply(lambda x: x.head(1))
display(top_cell_type_per_group.to_dict())

# 3dpi
top_cell_type_per_group_3 = day3_adata_gene_rank.groupby('group')['name'].apply(lambda x: x.head(1))
display(top_cell_type_per_group.to_dict())

sc.pl.umap(score_mock, color=["Ciliated cells","leiden_res_"], cmap="coolwarm")
sc.pl.umap(score_1, color=["Ciliated cells","leiden_res_"], cmap="coolwarm")
sc.pl.umap(score_2, color=["Ciliated cells","leiden_res_"], cmap="coolwarm")
sc.pl.umap(score_3, color=["Ciliated cells","leiden_res_"], cmap="coolwarm")    # Use "Clara cells" instead, since it appears at 3dpi

dict_ann_mock = mock_adata_gene_rank[mock_adata_gene_rank["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()
dict_ann_mock

dict_ann_1 = day1_adata_gene_rank[day1_adata_gene_rank["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()
dict_ann_1

dict_ann_2 = day2_adata_gene_rank[day2_adata_gene_rank["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()
dict_ann_2

dict_ann_3 = day3_adata_gene_rank[day3_adata_gene_rank["stat"] > 0].groupby("group").head(1).set_index("group")["name"].to_dict()
dict_ann_3

mock_adata.obs["leiden_res_"] = mock_adata.obs["leiden_res_"].cat.rename_categories(dict_ann_mock)
day1_adata.obs["leiden_res_"] = day1_adata.obs["leiden_res_"].cat.rename_categories(dict_ann_1)
day2_adata.obs["leiden_res_"] = day2_adata.obs["leiden_res_"].cat.rename_categories(dict_ann_2)
day3_adata.obs["leiden_res_"] = day3_adata.obs["leiden_res_"].cat.rename_categories(dict_ann_3)

sc.pl.umap(
    adata=mock_adata,
    color=[ "leiden_res_"],
    ncols=1,
    title='mock'
)

sc.pl.umap(
    adata=day1_adata,
    color=[ "leiden_res_"],
    ncols=1,
    title='1dpi'
)

sc.pl.umap(
    adata=day2_adata,
    color=[ "leiden_res_"],
    ncols=1,
    title='2dpi'
)

sc.pl.umap(
    adata=day3_adata,
    color=[ "leiden_res_"],
    ncols=1,
    title='3dpi'
)
```

## Trajectory Inference
It helps to understand how cells transition from one type to another, for example, stem cell to mature cell.

First, build a graph (network) connecting cells that look alike, where each cell is a dot, and a line (edge) drawn between two similar cells.
```py
# Trajectory analysis
sc.tl.draw_graph(mock_adata)
sc.tl.draw_graph(day1_adata)
sc.tl.draw_graph(day2_adata)
sc.tl.draw_graph(day3_adata)

plt.rcParams["figure.figsize"] = (4,4)
sc.pl.draw_graph(mock_adata, color='leiden_res_', size = 16, title='mock')
sc.pl.draw_graph(day1_adata, color='leiden_res_', size = 16, title='1dpi')
sc.pl.draw_graph(day2_adata, color='leiden_res_', size = 16, title='2dpi')
sc.pl.draw_graph(day3_adata, color='leiden_res_', size = 16, title='3dpi')
```
Next, ABSTRACT the graph, where all the points clustering for one cell type are converted to one point, Similar to a blunt summary of everypoint

```py
sc.tl.paga(mock_adata, groups='leiden_res_')
sc.tl.paga(day1_adata, groups='leiden_res_')
sc.tl.paga(day2_adata, groups='leiden_res_')
sc.tl.paga(day3_adata, groups='leiden_res_')

sc.pl.paga(mock_adata, color=['leiden_res_'], title='mock')
sc.pl.paga(day1_adata, color=['leiden_res_'], title='1dpi')
sc.pl.paga(day2_adata, color=['leiden_res_'], title='2dpi')
sc.pl.paga(day3_adata, color=['leiden_res_'], title='3dpi')

sc.tl.draw_graph(mock_adata, init_pos='paga')
sc.tl.draw_graph(day1_adata, init_pos='paga')
sc.tl.draw_graph(day2_adata, init_pos='paga')
sc.tl.draw_graph(day3_adata, init_pos='paga')

sc.pl.draw_graph(mock_adata, color='leiden_res_', legend_loc='on data', size=8, title='mock')
sc.pl.draw_graph(day1_adata, color='leiden_res_', legend_loc='on data', size=8, title='1dpi')
sc.pl.draw_graph(day2_adata, color='leiden_res_', legend_loc='on data', size=8, title='2dpi')
sc.pl.draw_graph(day3_adata, color='leiden_res_', legend_loc='on data', size=8,title='3dpi')

plt.rcParams["figure.figsize"] = (5,4)
sc.pl.paga_compare(mock_adata, threshold=0.03, frameon=True, edges=True, size = 16, title='mock')
sc.pl.paga_compare(day1_adata, threshold=0.03, frameon=True, edges=True, size = 16, title='1dpi')
sc.pl.paga_compare(day2_adata, threshold=0.03, frameon=True, edges=True, size = 16, title='2dpi')
sc.pl.paga_compare(day3_adata, threshold=0.03, frameon=True, edges=True, size = 16, title='3dpi')
```

Now how do cells transition from one type to another, assuming we have a pluripotent progenitor cell.

```py
mock_adata.uns['iroot'] = np.flatnonzero(mock_adata.obs['leiden_res_']  == 'Pluripotent stem cells')[0]
sc.tl.dpt(mock_adata)

day1_adata.uns['iroot'] = np.flatnonzero(day1_adata.obs['leiden_res_']  == 'Pluripotent stem cells')[0]
sc.tl.dpt(day1_adata)

day2_adata.uns['iroot'] = np.flatnonzero(day2_adata.obs['leiden_res_']  == 'Pluripotent stem cells')[0]
sc.tl.dpt(day2_adata)

day3_adata.uns['iroot'] = np.flatnonzero(day3_adata.obs['leiden_res_']  == 'Pluripotent stem cells')[0]
sc.tl.dpt(day3_adata)

sc.pl.draw_graph(mock_adata, color=['dpt_pseudotime', 'leiden_res_'], legend_loc='on data', size = 24, title='mock')
sc.pl.draw_graph(day1_adata, color=['dpt_pseudotime', 'leiden_res_'], legend_loc='on data', size = 24, title='1dpi')
sc.pl.draw_graph(day2_adata, color=['dpt_pseudotime', 'leiden_res_'], legend_loc='on data', size = 24, title='2dpi')
sc.pl.draw_graph(day3_adata, color=['dpt_pseudotime', 'leiden_res_'], legend_loc='on data', size = 24, title='3dpi')
```

Saving the datasets

```py
mock_adata.write("mock_adata.h5", compression="gzip")
day1_adata.write("day1_adata.h5", compression="gzip")
day2_adata.write("day2_adata.h5", compression="gzip")
day3_adata.write("day3_adata.h5", compression="gzip")
```
## Biological Interpretation: Questions and Answers
### Q1. What cell types did you identify at the different stages of infection?

### Q2. Can you explain why these cell types correlate with COVID-19 infection?

### Q3. Is ACE2 a good marker for tracking COVID-19 infection rate (based on this dataset)?

### Q4. What is the difference between ENO2 and ACE2 as biomarkers in the two studies?

### Q5. Which cell cluster has the highest abundance of ACE2 expression after 3 dpi and what does that mean biologically (INTERPRET VISUALLY)?

## REFERENCES
Ravindra, N.G., Alfajaro, m.M., Gasque, V., Huston, N.C., Wan, H., Szigeti-Buck, K., Yasumoto, Y., Greaney, A.M., Habet, V., Chow, R.D., Chen, J.S., Wei, J., Filler, R.B., Wang, B., Wang, G., Niklason, L.E., Montgomery, R.R., Eisenbarth, S.C., Chen, S., Williams, A., Iwasaki, A., Horvath, T.L., Foxman, E.F., Pierce, R.W., Pyle, A.M., van Dijk, D., Wilen, C.B. (2021). Single-cell longitudinal analysis of SARS-CoV-2 infection in human airway epithelium identifies target cells, alterations in gene expression, and cell state changes. _PLoS Biol_ 19(3): e3001143. https://doi.org/10.1371/journal.pbio.3001143.

