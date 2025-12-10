# Reproducing a Core Single-Cell RNA-Seq Analysis Pipeline Using Scanpy

The project aims at reproducing a **complete single-cell RNA-sequencing (scRNA-seq)** data analysis pipeline using the Python package **Scanpy**, covering the following steps:

1. Installation of packages
2. Loading and preprocessing of data  
3. Quality control (QC)  
4. Normalization and feature selection  
5. Dimensionality reduction (PCA, UMAP)  
6. Clustering  
7. Cell type annotation and visualisations
8. Biological interpretation of results

## Installation
All necessary Python packages were installed as follows:
```py
!pip install scanpy
!pip install anndata
!pip3 install igraph
!pip install celltypist              # Decoupler was eventually used instead
!pip install decoupler
!pip install fa2-modified
```

## Loading Data
In this step, the single-cell expression matrix was loaded into an **AnnData** object (the core data structure in Scanpy), which contains:
- `adata.X`: the expression matrix (cells × genes)
- `adata.obs`: metadata for each cell
- `adata.var`: metadata for each gene

```py
# Import core single cell tools
import scanpy as sc
import anndata as ad

# Download dataset
!wget https://github.com/josoga2/sc/raw/refs/heads/main/bone_marrow.h5ad

# Load and view dataset
bone_marrow_adata = sc.read_h5ad('bone_marrow.h5ad')
print(bone_marrow_adata)

# Print the dimensions of our dataset
bone_marrow_adata.shape               # 14783 cells and 17374 genes

# View the first 5 rows describing the genes in the dataset
bone_marrow_adata.var.head()

# Also view the first 5 rows describing the cells (ID) in the dataset
bone_marrow_adata.obs.head()

# View the dataset in a proper dataframe format
bone_marrow_adata.to_df()
```

## Quality Control (QC)
Ideally, QC ensures we only keep high-quality cells and informative genes, and typical filters remove:
- Harmonize unique gene names (avoid gene duplications from old pipelines)
- Cells with too few genes (likely dead)
- Cells with too many genes (possible doublets)
- Genes expressed in very few cells (uninformative)

```py
# Following a useful step for older datasets, just in case of necessity.
bone_marrow_adata.var_names_make_unique()
bone_marrow_adata.obs_names_make_unique()

# Checking for possible contamination from dying cells, ribosomal transcripts or haemoglobin

# As a rule of thumb:
# Cells with a high proportion of mitochondrial reads (say >10–20%) are likely stressed, apoptotic, or poorly captured.
# Ribosomal transcripts are removed because they represent global transcriptional activity, not cell-type-specific biology.
# Instead of true cell populations, high HB signal often represents ambient RNA contamination from lysed red blood cells.

bone_marrow_adata.var['MT'] = bone_marrow_adata.var_names.str.startswith("MT-")
bone_marrow_adata.var['RIBO'] = bone_marrow_adata.var_names.str.startswith("RPS", "RPL")
bone_marrow_adata.var['HB'] = bone_marrow_adata.var_names.str.startswith("^HB[^(P)]")

# Now a quick check.
mt_genes = bone_marrow_adata.var[bone_marrow_adata.var['MT']]
mt_genes

ribo_genes = bone_marrow_adata.var[bone_marrow_adata.var['RIBO']]
ribo_genes

hb_genes = bone_marrow_adata.var[bone_marrow_adata.var['HB']]
hb_genes
```
The quick check above showed there were no mitochondrial, ribosomal, or haemoglobin genes. Few more QC steps to ascertain this outcome:
```py
# calculate the qc metrics

sc.pp.calculate_qc_metrics(
    bone_marrow_adata, qc_vars=["MT", 'RIBO', 'HB'], inplace=True, log1p=True
)

# Note that it is also included in the headers of obs
bone_marrow_adata.obs.head()

# And the gene list
bone_marrow_adata.var.head()

# Average number of genes that have at least one detected identifier in each cell i.e., the number of genes expressed in each cell:
sc.pl.violin(
    bone_marrow_adata,
    ["n_genes_by_counts"],
    jitter=0.4,
    multi_panel=False,
)

# And the total number of molecules (UMI) detected in a cell:
sc.pl.violin(
    bone_marrow_adata,
    ["total_counts"],
    jitter=0.4,
    multi_panel=False,
)

# Are any mitochondrial genes present?
sc.pl.violin(
    bone_marrow_adata,
    ["pct_counts_MT"],
    jitter=0.4,
    multi_panel=False,
)

# What about the ribosomal genes?
sc.pl.violin(
    bone_marrow_adata,
    ["pct_counts_RIBO"],
    jitter=0.4,
    multi_panel=False,
)

# And the haemoglobin genes?
sc.pl.violin(
    bone_marrow_adata,
    ["pct_counts_HB"],
    jitter=0.4,
    multi_panel=False,
)

# Visualise each of them in scatter plots.
sc.pl.scatter(bone_marrow_adata, "total_counts", "n_genes_by_counts", color="pct_counts_MT")
sc.pl.scatter(bone_marrow_adata, "total_counts", "n_genes_by_counts", color="pct_counts_RIBO")
sc.pl.scatter(bone_marrow_adata, "total_counts", "n_genes_by_counts", color="pct_counts_HB")
```
All the steps above confirmed there no MT, RIBO, and HB genes in the dataset.
```py
sc.pp.calculate_qc_metrics(bone_marrow_adata, inplace=True, log1p=True)

# Reconfirming the dataset
bone_marrow_adata.obs.head()
bone_marrow_adata.var.head()

# Playing around the dataset with violin plots and scatterplots
sc.pl.violin(bone_marrow_adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4)
sc.pl.violin(bone_marrow_adata, ['n_genes_by_counts'], jitter=0.4)
sc.pl.violin(bone_marrow_adata, ['total_counts'], jitter=0.4)
sc.pl.scatter(bone_marrow_adata, x='total_counts', y='n_genes_by_counts', color='n_genes_by_counts')
sc.pl.scatter(bone_marrow_adata, x='total_counts', y='n_genes_by_counts', color='total_counts')
```
There was no need for the filtering step based on the QC outcomes.

### Doublet detection
```py
sc.pp.scrublet(bone_marrow_adata)
```

## Normalization
Since normalization adjusts for sequencing depth differences between cells, counts in the dataset were scaled so that each cell can the same total expression level.
```py
# First, save a copy of the data as follows:
bone_marrow_adata.layers["counts"] = bone_marrow_adata.X.copy()

# Normalizing to median total counts
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
### Q1. What cell types did you identify?
The following cell types were identified, based on each cluster annotation:

- Cluster 0: Neutrophils
- Cluster 1: Gamma delta T cells
- Cluster 2: T memory cells
- Cluster 3: NK cells
- Cluster 4: Nuocytes
- Cluster 5: B cells naive
- Cluster 6: Platelets
- Cluster 7: Plasma cells
- Cluster 8: Monocytes

### Q2. Explain the biological role of each cell type

**1. Neutrophils**: They are produced and stored in large pools within the bone marrow, from where they are rapidly mobilised. In peripheral tissues, they function as first-line innate defenders through phagocytosis, release of reactive oxygen species (ROS), and formation of neutrophil extracellular traps (NETs).

**2. Gamma Delta (γδ) T Cells**: They arise from bone-marrow–derived lymphoid progenitors but develop in the thymus. In peripheral immunity, γδ T cells provide rapid responses at epithelial and mucosal surfaces, recognising non-peptide antigens in an MHC-independent manner, and contribute early cytokines during infection, tumour or tissue stress.

**3. T Memory Cells**: They are subsets of T cells that home to the bone marrow for long-term maintenance after antigen-driven activation of naïve T cells. In the periphery, they provide accelerated and enhanced immune responses upon re-exposure to previously encountered antigens, a major feature of adaptive immunity.

**4. NK Cells (Natural Killer Cells)**: They originate from bone-marrow progenitors and circulate as cytotoxic innate lymphocytes. In peripheral immunity, they eliminate virus-infected and transformed (tumour) cells and produce IFN-γ to modulate immune responses.

**5. Nuocytes or Group 2 Innate Lymphoid Cells (ILC2s)**: They develop from innate lymphoid progenitors in the bone marrow. In peripheral tissues, they mediate type-2 immunity by producing IL-5 and IL-13, supporting anti-helminth defense, allergic inflammation, and tissue repair.

**6. Naïve B Cells**: Naïve B cells complete their development and B-cell receptor rearrangement in the bone marrow. After maturation, they enter peripheral circulation and localize to secondary lymphoid organs, where they survey for cognate antigen. Upon activation, naïve B cells differentiate into antibody-secreting plasma cells or memory B cells, initiating humoral immune responses.

**7. Platelets**: Platelets are generated by megakaryocytes in the bone marrow. In peripheral immunity, they function not only in haemostasis but also as immune mediators by releasing inflammatory factors and interacting with leukocytes. Their immunomodulatory abilities link coagulation with innate immunity.

**8. Plasma Cells**: Plasma cells represent the terminal differentiation stage of activated B cells. Long-lived plasma cells often home to the bone marrow, where they persist and continuously secrete high-affinity antibodies for long-term serological memory. Short-lived plasma cells in peripheral tissues produce immediate, high-titer antibody responses during acute infection. They are the primary mediators of humoral immunity.

**9. Monocytes**: They are produced from myeloid progenitors in the bone marrow and released into circulation. In peripheral immunity, they migrate into tissues, differentiate into macrophages or dendritic cells, and perform phagocytosis and antigen presentation.

### Q3. Is the tissue source really bone marrow? Justify your answer

It looks like the tissue source is NOT bone marrow for the following reasons:

**Missing bone marrow-specific progenitor populations**: The identified cell types in the clusters include (mature) neutrophils, monocytes, B and T lineages, NKs and platelets but no obvious haematopoietic progenitor or erythroid clusters that are bone marrow-specific.

**Typical peripheral blood frequency distribution of mature immune cells**: The dotplot, for example, shows many mature lymphoid (T memory, gamma delta T cells, NK cells, naive B cells, plasma cells) and myeloid (neutrophils, monocytes) populations at sizes and frequencies consistent with a peripheral blood sample. Bone marrow usually has a different composition, which includes higher proportions of immature myeloid precursors and erythroid cells, and often different relative abundances of mature circulating lymphocytes.

**Presence of platelets with other circulating leukocytes**:
Platelets were identified alongside mature leukocyte types, a pattern expected in whole blood or PBMC preps contaminated with platelets. In bone marrow, platelets are produced by megakaryocytes, hence, one would expect to see megakaryocyte lineage signatures rather than free platelets dominating a cluster.

In conclusion, it is more likely the tissue is a peripheral blood or whole blood sample (or PBMCs with platelet contamination), and not bone marrow, except in any of the following scenarios where bone marrow could masquerade as peripheral blood:

(a) the assay lacked progenitor markers,

(b) the experiment intentionally depleted progenitors, or

(c) there was an annotation error.

### Q4. Based on the relative abundance of cell types, is the patient healthy or infected?
The identified cell-type proportions are most consistent with an active infection or inflammatory state (most likely an acute innate response, classically bacterial) rather than a healthy steady state, due to the following reasons:

The results show that neutrophils have the largest relative abundance (over 90% fraction of cells in a group) and the highest mean or median expression in group, compared to other identified cell types. Such outcome is contrary to what one would expect for a healthy resting peripheral blood profile. An increased neutrophil fraction is therefore strong evidence for an ongoing innate inflammatory response, particularly in an acute bacterial infection, where neutrophilia (increase in the frequency and number of neutrophils) is the dominant signal.

Furthermore, the number of monocytes also seems to be elevated relative to a healthy baseline. Monocytosis (increase in the number of circulating monocytes) commonly accompanies acute infection and inflammation, and can reflect mobilisation of inflammatory monocyte subsets that move to tissues and differentiate to macrophages or dendritic cells (DCs). Hence, neutrophilia and monocytosis strengthens the interpretation of an active inflammatory or infectious process rather than an isolated lymphoid perturbation.

Moreover, the results show a non-trivial NK cell activation states, suggesting a non-antiviral response. Numerically increased NK cells with transcriptional evidence of activation, such as cytotoxic or IFN-γ signatures, would indicate innate antiviral or early anti-infectious activity.

Finally, despite a mix of adaptive immunity (indicated by the presence of lymphoid cells), innate myeloid cells dominated the immune response. Such phenomenon is commonly seen in acute severe infection and systemic inflammation where innate cells expand and lymphocytes marginate, deplete or are redistributed into tissues.
