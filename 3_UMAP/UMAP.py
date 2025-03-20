import scanpy as sc
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt


#for some reason running it as .py is way faster, we think because the console output
#is slowing down jupyter a ton
def perform_UMAP(): 
    
    print("Running... Takes a couple minutes")
    
    #reproducibility
    np.random.seed(42)

    data = sc.read("2_PCA/PCA_data.h5ad")

    #can tune PCs
    sc.pp.neighbors(data, n_pcs=50)

    sc.tl.umap(data, min_dist=0.3)

    #tune resolution based on diff gene expression
    resolution = 0.3
    sc.tl.leiden(data, flavor='igraph', n_iterations=-1, resolution=resolution) 
            
    sc.pl.umap(data, color='leiden', title='Leiden Clustering Results', show=False) #plot again with colours
    
    name = "3_UMAP/UMAP_Leiden_graph_" + str(resolution) + "_res.png"
    plt.savefig(name, dpi=300, bbox_inches="tight")  #save it
    
    name = "3_UMAP/UMAP_data_" + str(resolution) + "_res.h5ad"
    data.write(name)
    
    print(data.obs['leiden'].value_counts())
    
    return
    
if __name__ == "__main__":
    perform_UMAP()
    
