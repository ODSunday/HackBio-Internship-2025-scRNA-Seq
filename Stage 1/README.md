### Task 1: Essay
## Single-Cell Genomics: A Stepping Stone for Future Immunology Discoveries
The immune system is one of the most dynamic and adaptable systems in the human body, creating a complex network that connects nearly every tissue. Characterising the unique interactions of immune cells is crucial to enhancing therapeutic effectiveness. Since the 19th century, studies pioneered by Ilya Mechnikov have categorised immune cells into distinct types based on their morphology and cellular functions. Technological advancements have accelerated immunological identification studies. However, the belief that a small number of molecular determinants sufficiently describe the diversity of immune cells has been disproven, as immune cell identity largely depends on their tissue and environmental context. 

Single-cell sequencing, therefore, solves the challenge associated with the heterogeneity of immune cells, as it can characterise the various cell types and states involved in immune cell functions. Single-cell sequencing facilitates the modelling of gene expression profiles from a bulk of cells, thereby overcoming the limitations of bulk sequencing, including challenges in reproducibility, sensitivity, automation, microfluidic technologies, and noise. 

Key challenges are limiting the full potential of single-cell genomics in immunology. These include analytical complexity occasioned by large datasets that encompass thousands of cells having diverse transcriptional profiles, and other computational hassles such as accurate classification of cell types, differentiating subtle functional states, and filtering technical noise. 

Importantly, single-cell approaches now integrate with flow cytometry, CyTOF, imaging, lineage tracing, spatial transcriptomics, and CRISPR-based perturbation assays. These multimodal strategies reconnect transcriptional profiles with protein expression, spatial location, and clonal history, promising to restore contextual understanding while overcoming limitations of dissociated single-cell RNA-sequencing. Collectively, these advances reveal the full diversity, organisation, and functional states of the immune system across health and disease. Hence, single-cell genomics is an indispensable tool for decoding immune diversity, playing a major role in biomarker discovery and therapeutic innovation. 

**Keywords**: Single-cell genomics, immune diversity, immunology, therapeutic innovation, and biomarker.

### Task 2: Plotting
## Part A. Gene Expression (Heatmap & Volcano Plot)
a. Normalized counts for HBR vs UHR samples

The normalized gene expression dataset was used to plot a clustered heatmap of the top differentially expressed genes (DEGs) between HBR and UHR samples.
- Labels: genes and samples.
- Colour gradient (Blues) indicates expression levels.
  
Dataset: https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv
```py
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Import dataset
data_source = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv"
df = pd.read_csv(data_source, index_col=0)

# Inspect the data
print(df.shape)
print(df.head())

# Create a clustered heatmap
sns.set(style="white", font_scale=0.9)
g = sns.clustermap(
    df,
    cmap="Blues",        # colour gradient
    linewidths=0.7,
    linecolor='black',
    figsize=(4, 4),
    xticklabels=True,
    yticklabels=True
)

plt.show()               # Shows plot
```

#### Resulting plot

a.

<img width="395" height="395" alt="figure 1a" src="https://github.com/user-attachments/assets/dc48c78a-19fb-4dc7-bc4e-3ed0e0997bae" />

###### Figure 1a: Heatmap of the top DEGs between HBR and UHR samples

b. Differential expression results (chromosome 22)
The plot of log2FoldChange vs log10(Padj) from the DEG results.
- Colour points by significance:
  - Upregulated: `green`
  - Downregulated: `orange`
  - Not significant: `grey`
- Dashed vertical lines at log2FoldChange = ±1.
  
Dataset: https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv
```py
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Import dataset
data_source = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv"
df = pd.read_csv(data_source, index_col=0)

# Inspect the data
print(df.shape)
print(df.head())
print(df.tail())

# Compute -log10 adjusted p-value
df["-log10PAdj"] = -np.log10(df["PAdj"])

# Define significance categories
def gene_category(row):
    if (row["PAdj"] < 0.05) and (row["log2FoldChange"] > 1):
        return "Upregulated"
    elif (row["PAdj"] < 0.05) and (row["log2FoldChange"] < -1):
        return "Downregulated"
    else:
        return "Not significant"

df["Significance"] = df.apply(gene_category, axis=1)

# Set colour palette
palette = {
    "Upregulated": "green",
    "Downregulated": "orange",
    "Not significant": "grey"
}

# Create the Volcano Plot
plt.figure(figsize=(8, 6))
sns.scatterplot(
    data=df,
    x="log2FoldChange",
    y="-log10PAdj",
    hue="Significance",
    palette=palette,
    alpha=0.7,
    edgecolor=None
)

# Add significance threshold lines
plt.axvline(x=1, color='black', linestyle='--', linewidth=1)
plt.axvline(x=-1, color='black', linestyle='--', linewidth=1)
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', linewidth=1)

# Label and style
plt.title("Volcano Plot of log2FoldChange vs -log10PAdj", fontsize=12)
plt.xlabel("log2FoldChange", fontsize=10)
plt.ylabel("-log10PAdj", fontsize=10)
plt.legend(title="Gene Regulation", loc="upper right")
sns.despine()

plt.tight_layout()
plt.show()
```

#### Resulting plot

b.

<img width="787" height="587" alt="figure 1b" src="https://github.com/user-attachments/assets/28353c76-2f75-43a7-b36a-019cc5bd5fdd" />

###### Figure 1b: Volcano Plot of log2FoldChange vs -log10PAdj

## Part B. Breast Cancer Data Exploration
c. Scatter Plot (radius vs texture)

The scatter plot of `texture_mean` vs `radius_mean`: Colour points by `diagnosis` (M = malignant, B = benign).

Dataset: https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv

```py
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Import dataset
data_source = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
df = pd.read_csv(data_source, index_col=0)

# Inspect the data
print(df.shape)
print(df.head())
print(df.tail())

# Check for NaN values in the dataset
print(df.isnull().sum())

# Relevant columns do not contain missing values. Proceed to the next steps.

# Create the scatter plot
plt.figure(figsize=(6, 5))
sns.scatterplot(data=df, x='radius_mean', y='texture_mean', hue='diagnosis', palette={'M': 'blue', 'B': 'orange'})

# Add titles and labels
plt.title('Scatter Plot of Texture Mean vs Radius Mean')
plt.xlabel('radius_mean')
plt.ylabel('texture_mean')
plt.legend(title='diagnosis', loc='upper right')

# Show the plot
plt.show()
```

#### Resulting plot

c.

<img width="536" height="474" alt="figure 1c" src="https://github.com/user-attachments/assets/dd109aae-d33b-40d4-a5ae-856d0606faa1" />

###### Figure 1c: scatter plot of texture_mean vs radius_mean

d. Correlation Heatmap
Compute the correlation matrix of six key features:

- `radius_mean`, `texture_mean`, `perimeter_mean`, `area_mean`, `smoothness_mean`, `compactness_mean`.

Plot as a heatmap with correlation values annotated.

Dataset: https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv
```py
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Import dataset
data_source = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
df = pd.read_csv(data_source, index_col=0)

# Inspect the data
print(df.shape)
print(df.head())
print(df.tail())

# Select the six key features
key_features = ['radius_mean', 'texture_mean', 'perimeter_mean', 'area_mean', 'smoothness_mean', 'compactness_mean']
df_selected = df[key_features]

# Compute the correlation matrix
correlation_matrix = df_selected.corr()

# Set up the matplotlib figure
plt.figure(figsize=(8, 6))

# Create a heatmap with annotations
sns.heatmap(correlation_matrix, annot=True, fmt=".1f", cmap='Blues', square=True, cbar_kws={"shrink": .8})

# Add titles and labels
plt.title('Correlation Heatmap of Key Features', fontsize=12)
plt.xticks()
plt.yticks()

# Show the plot
plt.tight_layout()  # Adjusts layout to accommodate the labels
plt.show()
```

#### Resulting plot

d.

<img width="648" height="590" alt="figure 1d" src="https://github.com/user-attachments/assets/43d1bf0e-646b-4a23-9107-75a17dee0797" />

###### Figure 1d: Correlation heatmap of key features

e. Scatter Plot (smoothness vs compactness)

- Plot `compactness_mean` vs `smoothness_mean` coloured by `diagnosis`.
- Include gridlines and clear axis labels.

Dataset: https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv
```py
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Import dataset
data_source = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
df = pd.read_csv(data_source, index_col=0)

# Inspect the data
print(df.shape)
print(df.head())
print(df.tail()

# Create the scatter plot
plt.figure(figsize=(6, 5))
palette = {'M': 'blue', 'B': 'orange'}      # Defines colour palette
sns.scatterplot(data=df, x='smoothness_mean', y='compactness_mean',
                hue='diagnosis', palette=palette, alpha=0.7)

plt.grid(True)                              # Adds gridlines

# Add titles and labels
plt.title('Scatter Plot of Smoothness vs Compactness', fontsize=12)
plt.xlabel('smoothness_mean', fontsize=10)
plt.ylabel('compactness_mean', fontsize=10)
plt.legend(title='diagnosis', loc='upper left')

plt.tight_layout()                          # Adjusts layout to accommodate the labels
plt.show()                                  # Shows the plot
```

#### Resulting plot

e.

<img width="590" height="490" alt="figure 1e_2" src="https://github.com/user-attachments/assets/a71d3b93-02c8-40fe-99f1-a34b0f80f52b" />

###### Figure 1e. Scatter plot of compactness_mean vs smoothness_mean

f. Density Plot (area distribution)

- Plot kernel density estimates (KDE) of `area_mean` for both M and B diagnoses on the same axis.
- Add legend and labeled axes.

Dataset: https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv
```py
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Import dataset
data_source = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"
df = pd.read_csv(data_source, index_col=0)

# Inspect the data
print(df.shape)
print(df.head())
print(df.tail()

# Create the density plot using Seaborn
plt.figure(figsize=(8, 4))
sns.kdeplot(data=df[df['diagnosis'] == 'M']['area_mean'],
             color='blue', label='Malignant (M)', fill=True, alpha=0.5)
sns.kdeplot(data=df[df['diagnosis'] == 'B']['area_mean'],
             color='orange', label='Benign (B)', fill=True, alpha=0.5)

# Add titles and labels
plt.title('Kernel Density Estimate of Area Mean by Diagnosis', fontsize=12)
plt.xlabel('area_mean', fontsize=10)
plt.ylabel('Density', fontsize=10)

plt.legend()            # Adds legend

plt.tight_layout()      # Adjusts layout to accommodate the labels
plt.show()              # Shows the plot
```

#### Resulting plot

f.

<img width="789" height="390" alt="figure 1f" src="https://github.com/user-attachments/assets/c141a678-a97f-47b0-a0ee-079a8ee7583d" />

###### Figure 1f: Density plot of area distribution

### Task 3: DNA to protein translation
## Python function for translating DNA to protein

```py
# Starts at the first 'ATG' (start codon) and stops at the first stop codon, representing normal cellular process.

def translate_dna_to_protein(dna_seq):
    codon_table = {
        'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
        'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
        'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
        'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
        'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
        'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
        'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
        'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
        'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
        'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
        'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
        'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
        'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
        'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
        'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*',
        'TGC':'C','TGT':'C','TGA':'*','TGG':'W'
    }
    # Stop codons (TAA, TAG, TGA) are represented by "*".
    # Invalid or incomplete codons are translated as "X" to indicate an unknown amino acid.

    # Ensure uppercase and remove any space in the sequence string
    dna_seq = dna_seq.upper().replace(" ", "")

    # Find the first start codon (ATG)
    start_index = dna_seq.find("ATG")
    if start_index == -1:
        return "No start codon found."

    protein = ""

    # Translate starting from the first ATG.
    for i in range(start_index, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        amino_acid = codon_table.get(codon, 'X')
        if amino_acid == '*':  # stop codon
            break
        protein += amino_acid

    return protein


# DNA sequence (from the human AKT1 gene)
dna_sequence = (
    "ATGAGCGACGTGGCTATTGTGAAGGAGGGTTGGCTGGGCCCGAGTGGAAGGACAAGGGGCTGGAGGAGGAG"
    "CAGCAGGAGGCCATGCTTGGGGAAGGAGGAGGAGGAGGGCTGGAGCAGGAGGAGGAAGAGGCTGGAGGAAGGAG"
)

protein = translate_dna_to_protein(dna_sequence)
print(protein)
```

**Result**: `MSDVAIVKEGWLGPSGRTRGWRRSSRRPCLGKEEEEGWSRRRKRLEEG`

A quick `BLAST` of `MSDVAIVKEGWLGPSGRTRGWRRSSRRPCLGKEEEEGWSRRRKRLEEG` in `UniProtKB` identified the protein sequence as `non-specific serine/threonine protein kinase` (AKT1 gene)

### Task 4: Calculation of hamming distance
## Python function for calculating the hamming distance between slack username and twitter/X handle.

Hamming distance = number of differing positions
```py
def hamming_distance(handle1, handle2):                       # Defines the function
    distance = sum(a != b for a, b in zip(handle1, handle2))  # Compares characters in overlapping range
    return distance

slack = "Sunday"
twitter = "SundOD"                                            # Synthesised Twitter handle.

print("Hamming distance:", hamming_distance(slack, twitter))
```

**Result**: 2
