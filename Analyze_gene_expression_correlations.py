# import libraries
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import gseapy as gp
import numpy as np
from gseapy import barplot
from Bio import SeqIO
from scipy import stats
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from bioservices import KEGG

# Setup path for the genes
root_dir = Path(r'E:\Coding\Python\Biopython\biopython projects\Analyze gene expression correlations\data')
gene_folders = ['AURKA', 'BUB1', 'CCNB2', 'CDK1', 'FBXO5', 'KIF4A', 'MCM10', 'PLK1', 'TOP2A', 'TPX2']

all_data = []

# Loop through all the genes
for folder in gene_folders:
    fna_path = root_dir / folder / 'ncbi_dataset' / 'data' / 'gene.fna'

    if fna_path.exists():
        for record in SeqIO.parse(fna_path, "fasta"):
            seq_str = str(record.seq)
            all_data.append({
                'gene': folder.split('_')[0],
                'length': len(seq_str),
                'gc_content': (seq_str.count('G') + seq_str.count('C')) / len(seq_str) * 100
            })
    else:
            print(f"file fna_path not identified")

# Creating DataFrame and Print it out
df = pd.DataFrame(all_data)

if not df.empty:
    print("Data processed successfully!")
    print(df.head())

# Generate the BoxPlot
    plt.figure(figsize=(8,5))
    sns.boxplot(x='gene', y='gc_content', data=df)
    plt.title('GC Content Distribution by Gene')

# Group genes for ANOVA
    groups = [group['gc_content'].values for name, group in df.groupby('gene')]
    f_stat, p_val = stats.f_oneway(*groups)

    print(f"\n--- statistical validation---")
    print(f"ANOVA p-value: {p_val:.4e}")

    if p_val < 0.05:
         print("Significant difference in GC content across the selected genes")
    else:
         print("No significant difference in the GC content across the selected genes.")
    plt.show()
else:
    print("Error: no data retrieved, check folder structure")

# data preparation
features = df[['length', 'gc_content']]
x = StandardScaler().fit_transform(features)

# PCA
pca = PCA(n_components=2)
principalComponents = pca.fit_transform(x)
pca_df = pd.DataFrame(data=principalComponents, columns=['PCA1', 'PCA2'])
pca_df['gene'] = df['gene']

# Visualize
plt.figure(figsize=(8,6))
sns.scatterplot(x='PCA1', y='PCA2', hue='gene', data=pca_df, s=100)
plt.title('PCA Analysis of Gene Features')
plt.show()

# craeting local mapping dictionary based on known gene pathway links
gene_to_pathway = {
     'CCNB2': 'hsa04110', #cell cycle
     'KIF4A': 'hsa04114', #oocyte meiosis
     'MCM10': 'hsa04110', #cell cycle
     'CDK1': 'hsa04110', #cell cycle
     'AURKA': 'hsa04110', #cell cycle
     'TOP2A': 'hsa04110', #cell cycle
     'BUB1': 'hsa04110', #cell cycle
     'PLK1': 'hsa04110', #cell cycle
     'TPX2': 'hsa04110', #cell cycle
}

# Apply this to your DataFrame
df['pathway'] = df['gene'].map(gene_to_pathway)
df['pathway_name'] = df['pathway'].replace({
    'hsa04110': 'Cell Cycle',
    'hsa04114': 'Oocyte Meiosis'
})

print(df[['gene', 'pathway_name']])

# Gene Ontology Enrichment Analysis

names= gp.get_library_name()
matches = [name for name in names if 'GO_Biological_Process' in name]
print("Valid Libraby names found:", matches)

valid_libs = gp.get_library_name()
target_lib = 'GO_Biological_Process_2023'

if target_lib in valid_libs:
     enr = gp.enrichr(gene_list=['AURKA', 'BUB1', 'CCNB2', 'CDK1', 'FBXO5', 'KIF4A', 'MCM10', 'PLK1', 'TOP2A', 'TPX2'],
                     gene_sets=[target_lib],
                     organism = 'human',
                     outdir='result')
     print("Enrinchment Analysis DONE!")
     print(enr.results.head())
else:
     print(f"Library {target_lib} not found. please choose from: {valid_libs[:5]}")

sig_results = enr.results[enr.results['Adjusted P-value'] < 0.05].copy()
test_dir = Path('result')
test_dir.mkdir(exist_ok=True)


if not sig_results.empty:
    top_10 = sig_results.sort_values('Adjusted P-value').head(10)

    plt.figure(figsize=(10, 6))
    plt.barh(top_10['Term'], -np.log10(top_10['Adjusted P-value']), color='skyblue')

    plt.xlabel('-log10(Adjusted P-value)')
    plt.ylabel('Biological Process')
    plt.title('Top 10 Enriched GO Terms')
    plt.gca().invert_yaxis() 

    save_path = test_dir / 'enrichment_analysis.png'
    plt.savefig(save_path, bbox_inches='tight')
    print(f"Figure successfully saved to: {save_path}")
else:
     print("no significant pathways found.")