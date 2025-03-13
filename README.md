# TRIM-IT

TRIM-IT is an unsupervised variable selection tool designed to retain overall data information for exploratory analysis. In Coletti et al. (doi: …. ), TRIM-IT is applied to GBM transcriptomics data to reduce data dimensionality and stratify patients into reliable clusters. The methodology consists of four sequential steps, each implemented in separate scripts:
 
1-MOMIP-GBM:
* Applies the Max-Out Min-In Problem (MOMIP) to GBM data for an initial unsupervised variable selection.
* Results are saved in the MOMIP-GBM-results.RData file, located in the Results folder.
* Input required: RNAseq-protein_informed-2021WHO.RData from the Data folder. 

2-K-means-Clustering:
* Determines the optimal number of clusters (k) and applies K-means clustering to group GBM patients into k subgroups considering:
    1. The complete GBM dataset.
    2. The GBM dataset reduced by MOMIP.
* The two clustering results are compared and evaluated by silhouette and Calinski-Harabsz scores.
* Clusters are visualized using UMAP (Uniform Manifold Approximation and Projection), a dimensionality reduction technique that preserves the local and global structure of the data.
* Results are saved in the Kmeans-GBM-results.RData file, stored in the Results folder.
* Input required: RNAseq-protein_informed-2021WHO.RData from the Data folder; MOMIP-GBM-results.RData from the Results folder. 

3-Analysis-Clinical-Features:
* Comprises two sections:
    * Correlation analysis between GBM clusters and clinical features (e.g., histology, sex).
    * Survival analysis using Kaplan-Meier curves, with pairwise log-rank tests and Cox regression analysis for comparing survival curves across clusters. These methods are commonly used to estimate survival probabilities and compare survival distributions between patient subgroups.
* Input required: Glioma-clinic-TCGA.RData from the Data folder; Kmeans-GBM-results.RData from the Results folder.
  
4-DGE_Analysis:
* Performs differential gene expression analysis using the edgeR R package.
* Generates violin plots for differentially expressed genes, providing insights into the distribution and expression levels of selected genes across clusters.
* Input required: raw-counts-GBM.csv, raw-counts-LGG.csv, and RNAseq-protein_informed-2021WHO.RData from the Data folder; MOMIP-GBM-results.RData and Kmeans-GBM-results.RData from the Results folder.

To reproduce the analysis, please download all the content of this folder and save it into the same location.

