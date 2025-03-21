We analysed data from a single cell sequencing dataset, from cells in the malpighian tubules of a fly.

The structure of our repository is as follows (each file should be run in order):

-0_data : our initial dataset

-1_preprocessing: where we preprocess our data

-2_PCA: we perform PCA

-3_UMAP: we perform a UMAP (not an ipynb for technical reasons)

-4_diff_gene_analysis: we analyse our data in more detail in multiple notbooks (run dotplot.ipynb first, then marker_UMAP, then top_genes_to_txt)

-5_conclusion our final conclusions

Summary of the analysis (each file has a short paragraph detailing our thought process, that we're putting here for ease of use):

### The first step of the project is the preprocessing of the dataset, we want to remove aberrant cells/genes, so that we have a 'cleaner' output, thus improving the quality of our downstream analysis.

For quality control purposes, we display a histogram of the number of *unique* transcripts per cell, to see if some cell contain no unique transcripts, or if some are very unique, we want to get an idea of the heterogeneity of our dataset. Then we display a violon plot of the number of genes per cell, as well as the total transcripts per cell, to check for aberrant values, and to guide our filtering.

In the histogram and violin plots, we see a small outlier exhibiting significantly more genes, transcripts than the others. So we chose to remove it out of an abundance of caution. To keep only high quality reads, we filter out cells expressing less than 500 genes, and filter genes expressed in fewer than 100 cells. We also normalize in cpm and log + 1 to account for differences in total number of transcripts per cell.

After the processing, we see that the data is more homogenous. We lost 800 cells out of 7000 initially (~10%), and 4400 genes out of 12400 initially (~35%). Losing 10% of cells is fine, and losing 35% of genes is a little high but, that is fine since our threshold isn't extremely high (gene expressed in at least 100 cells).

### The second step in our pipeline is PCA, to reduce the dimensionality of our data, making it easier to run a UMAP in the next step and allowing us to check that the PCA looks normal

We see that the first PC by itself explains 40% of the variance, but its followed by sharp diminishing returns, with 50 PCs, we get an explained variance of ~65%, which if fine for single cell analysis. The high dimensionality of our dataset makes it hard to compress down into just 50 PCs in our case. `<br>`

We note the first two PCs cannot capture the data well by themeselves. We also notice that the points are smeared along the PC1 axis, reflecting its high.er explanatory power `<br>`

This justifies our choice of only using 50PCs for the downstream analysis since we have greatly diminishing returns on the explained variance.

### Next we run an UMAP

The UMAP allows us to visualize the highly dimentionnal data in 2D. The leiden algorithm allowws us to identify separate clusters and colour them accordingly.

Following this procedure, we obtain 11 distinct clusters. However, some contains very few cells (<100).

A set seed was used for reproducibility. We used a bigger than normal minimal distance forthe UMAP (0.3 vs. 0.1) so we could see the clusters better. We decided on using a resolution of 0.3 after trying multiple resolutions based on how the UMAP looked, and based on our downstream analysis. We still get clusters on the UMAP that look like they might be able to be combined, we will further explore this in the later analyses.

The UMAP is saved as 'UMAP_Leiden_graph_0.3_res.png

### Differential gene analysis

*Dotplot:*

The dotplot shows us, by columns, for each group, the top 5 differentially expressed genes, and their expression level across all groups (in the rows). The size of the dot shows the fraction of cells expressing that gene, and the colour shows the mean expression level in the group. `<br>`

Looking at group 0, the top genes between 0 and 8 are very similar, suggesting that we may be able to group them together. Looking at group 2 and 3, they also share similarities, suggesting that we may be able to group them together also. The rest of them look good, indeed ignoring these clusters, we see a clear diagonal showing that we have well defined clusters with different marker genes, indicating different cell types in each cluster.

*UMAP + violin:*

Looking at the results of the differential gene expression mixed with violin plots, we see that once again (dotplot), groups 0 and 8 share the same top marker genes. Group 1 doesn't have markers that are truly exclusive, but they are more highly expressed in group 1. Group 2 does share marker genes with group 3, but it also has a unique marker, justifying that its a seperate cluster, so we should not combine it with group 2! The rest of the groups all have unique marker genes, justifying that they are indeed seperate clusters of cell types.

This is also a recap of the names of the top marker genes for each cluster.

*Top marker genes:*

https://www.flyrnai.org/tools/single_cell/web/

We used this website to lookup the cell types based on marker genes.

Cluster 0: Lower tubule PC 93% overlap `<br>`

Cluster 1: Principle Cell 99% overlap `<br>`

Cluster 2: Tubule bar shaped cell of initial segment 79% overlap `<br>`

Cluster 3: Stellate 84% overlap `<br>`

Cluster 4: Initial/transitional PC 67% overlap `<br>`

Cluster 5: Main segment PC 58% overlap `<br>`

Cluster 6: Adult Renal Stem Cell 81% overlap `<br>`

Cluster 7: unidentified 79% overlap `<br>`

Cluster 8: Lower tubule PC 61% overlap `<br>`

Cluster 9: Tracheal cell 57% overlap `<br>`

Cluster 10: Principle cell of lower ureter 34%; posterior midgut 39% overlap `<br>`

We see that most of our clusters have different defined cell types (except for 0 and 8, further evidence that they should be grouped). Cluster 7 shows an unidentified cell type, and cluster 9 surprisingly shows tracheal cells, which is strange since those kinds of cells shouldn't appear in the tubules.

### Conclusion, we put our findings together

We see 11 distinct cell types, that have been succesfully seperated, each with distinct marker genes, most of the cell types make sense for the malpighian tubules, except for the tracheal cells and the unidentified cluster. The UMAP is saved as UMAP_Leiden_graph_final.png
