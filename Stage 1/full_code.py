### Part A. Gene Expression (Heatmap & Volcano Plot)
## a. Normalized counts for HBR vs UHR samples

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


## b. Differential expression results (chromosome 22)

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


### Part B. Breast Cancer Data Exploration
## c. Scatter Plot (radius vs texture)

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

#d. Correlation Heatmap
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

## e. Scatter Plot (smoothness vs compactness)

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

## f. Density Plot (area distribution)

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

## Python function for translating DNA to protein

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
print(protein)          # Prints MSDVAIVKEGWLGPSGRTRGWRRSSRRPCLGKEEEEGWSRRRKRLEEG


## Python function for calculating the hamming distance between slack username and twitter/X handle

def hamming_distance(handle1, handle2):                       # Defines the function
    distance = sum(a != b for a, b in zip(handle1, handle2))  # Compares characters in overlapping range
    return distance

slack = "Sunday"
twitter = "SundOD"                                            # Synthesised Twitter handle.

print("Hamming distance:", hamming_distance(slack, twitter))  # Hamming distance: 2
