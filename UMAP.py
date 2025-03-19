import scanpy as sc
import anndata as ad
import numpy as np
import os 
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
#python -W ignore UMAP.py

def plot_qc_metrics(data):

    sns.histplot(data.obs["n_genes_by_counts"].unique(), bins=100)
    plt.title("Unique transcripts per cell")

    plt.tight_layout()
    plt.show()
    
    sc.pl.violin(data, ['n_genes_by_counts','total_counts'], multi_panel=True, show=False)
    plt.gcf().axes[0].set_title('Number of Genes by Counts')
    plt.gcf().axes[1].set_title('Total transcripts per cell')

    plt.tight_layout()
    plt.show()


#reproducibility
np.random.seed(42)

output_dir = "./Output"
file_path = "./data/malpighian_tubule.tsv"
os.makedirs(output_dir, exist_ok=True)

data = sc.read(file_path)

data = data.transpose() #need to transpose it to get it in the right way scanpy expects (cell, genes)

print(f"Number of cells: {data.n_obs}")
print(f"Number of genes: {data.n_vars}")

sc.pp.calculate_qc_metrics(data, inplace=True)
plot_qc_metrics(data)

sc.pp.filter_cells(data, min_genes=500)  #filter out cells with fewer than 500 genes
sc.pp.filter_genes(data, min_cells=100)  #filter out genes expressed in fewer than 10 cells
data = data[data.obs["n_genes_by_counts"] <= 5500].copy()

print(f"Number of cells after filtering: {data.n_obs}")
print(f"Number of genes after filtering: {data.n_vars}")

sc.pp.normalize_total(data.copy(), target_sum=1e6) #normalize in cpm
sc.pp.log1p(data) #log normalize

plot_qc_metrics(data) #We have one group of cell with high number of genes, look at where these cells are on the PCA to see if this is an artifact or not, are these just one custer

sc.tl.pca(data, svd_solver='arpack', n_comps = 50)

#elbow plot
plt.figure(figsize=(8, 5))
plt.plot(range(1, len(data.uns['pca']['variance_ratio']) + 1), data.uns['pca']['variance_ratio'], marker='o')
plt.xlabel('Principal Component')
plt.ylabel('Proportion of Variance Explained')
plt.title('Elbow Plot for PCA')
plt.grid(True)
plt.show()

#2D plot
sc.pl.pca(data, title='2D PCA')

#can tune PCs
sc.pp.neighbors(data, n_pcs=50)

sc.tl.umap(data, min_dist=0.3)

#tune resolution based on diff gene expression
sc.tl.leiden(data, flavor='igraph', n_iterations=-1, resolution=0.3) 
sc.pl.umap(data, color='leiden', title='Leiden Clustering Results')


sc.tl.rank_genes_groups(data, groupby="leiden", method="wilcoxon")
groups = data.obs["leiden"].unique()
sc.pl.rank_genes_groups(data)


top_genes=[]
for group in groups:
    # Get the ranked genes for each group
    df = sc.get.rank_genes_groups_df(data, group)
    top_genes_group = df["names"].head(3).tolist()  # Change 10 to the desired number of top genes

    top_genes.append(top_genes_group)

    df.to_csv(f"genes_{group}.csv")
top_genes = [gene for sublist in top_genes for gene in sublist]

#sc.pl.rank_genes_groups_dotplot(data, groupby="leiden", standard_scale="var", n_genes=5)

#sc.pl.heatmap(data, groupby='leiden')

output_dir = "./Output"
os.makedirs(output_dir, exist_ok=True)

# List of genes to visualize

for gene in top_genes:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # here we'll colour it by expression of the gene
    sc.pl.umap(data, color=gene, ax=axes[0], show=False, color_map="viridis")
    axes[0].set_title(f"UMAP - {gene} expression")

    # Violin plot of gene expression across Leiden clusters
    sc.pl.violin(data, keys=gene, groupby="leiden", ax=axes[1], show=False)
    axes[1].set_title(f"Violin Plot - {gene} expression")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{gene}_UMAP_Violin.png"), dpi=300)
    plt.close()

print(f"Plots saved in {output_dir}")

data.write("processed_data.h5ad")  

data = sc.read("processed_data.h5ad")
