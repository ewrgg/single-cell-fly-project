import scanpy as sc
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt


#Rnning it as .py is y faster, as the console output is slowing down jupyter due to 
#a warning caused by dependencies conflict in matplotlib version
def perform_UMAP(): 
    
    print("Running... Takes a couple minutes. If console fills with errors, thats normal.")
    
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
    
"""
The UMAP allows us to visualize the highly dimentionnal data in 2D. 
The leiden algorithm allowws us to identify separate clusters and colour them accordingly. 
Following this procedure, we obtain 11 distinct clusters. However, some contains very few cells (<100).

A set seed was used for reproducibility. We used a bigger than normal minimal distance forthe UMAP (0.3 vs. 0.1) 
so we could see the clusters better. We decided on using a resolution of 0.3 after trying multiple resolutions
based on how the UMAP looked, and based on our downstream analysis. We still get clusters on the UMAP that
look like they might be able to be combined, we will further explore this in the later analyses.

The UMAP is saved as 'UMAP_Leiden_graph_0.3_res.png
"""