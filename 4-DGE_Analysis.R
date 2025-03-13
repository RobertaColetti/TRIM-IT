### Step 4 of TRIM-IT: DIFFERENTIAL GENE EXPRESSION (DGE) ANALYSIS

# Library section
library(readr)
library(edgeR)
library(ggplot2) 
library(dplyr)

# Data loading and preparing.
#
# Load raw data for the DGE analysis
load("~/Data/rawCounts_GBM.RData") #here samples are on the columns
# Load MOMIP results to focus the analysis on the variables of interest
load("~/Results/MOMIP-GBM-results.RData")
variables=var_sel.2
# Load the reference dataset for graphical representation
load("~/Data/RNAseq-protein_informed-2021WHO.RData")

#Prepare data for the analysis
rawGBM_data=rawGBM_data[which(rownames(rawGBM_data) %in% variables),] #filter only the variables of interest
rawGBM_data=rawGBM_data[,order(colnames(rawGBM_data))]

# Load clustering results
load("~/Results/Kmeans-GBM-results.RData")
clusters=as.data.frame(res.km.MOMIP$cluster)
colnames(clusters)="cluster"
all(rownames(clusters)==colnames(rawGBM_data)) #it MUST be TRUE


## .............. DGE Analysis..............
#
set.seed(123)
y=estimateDisp(rawGBM_data,design = clusters) #common/trended/tagwise dispersions (necessary for next step)
design <- model.matrix(~factor(clusters$cluster))
#Perform quasi-likelihood F- test for DGE
fit <- glmQLFit(rawGBM_data, design,dispersion = y$trended.dispersion)

# Results
#C2 vs C1: 
set.seed(123)
qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1)

# C1 vs C3:
set.seed(123)
qlf.3vs1 <- glmQLFTest(fit, coef=3)
topTags(qlf.3vs1)

#C2 vs C3: 
set.seed(123)
qlf.2vs3 <- glmQLFTest(fit, contrast=c(0,-1,1))
topTags(qlf.2vs3)

#  .......... Construction of Violin Plot  ...........
# Prepare data for plotting (here normalized transcriptomics data are needed)
expression_data_norm <- as.data.frame(cbind(gene = colnames(GBM_RNA), t(GBM_RNA))) # For visualizing and comparing the gene expression among clusters, normalized data are preferable
expression_data_norm=expression_data_norm[,order(colnames(expression_data_norm))]
  sample_info <- as.data.frame(cbind(
      sample = colnames(rawGBM_data),
      Clusters = clusters$cluster
  ))
# Convert in long format for the plot
expression_long <- expression_data_norm %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  left_join(sample_info, by = "sample")
# Filter genes of interest identified from DGE analysis
C2vsC1=topTags(qlf.2vs1, n = Inf)$table
C1vsC3=topTags(qlf.3vs1, n = Inf)$table
C2vsC3=topTags(qlf.2vs3, n = Inf)$table
genes_of_interest = c(rownames(C2vsC1)[C2vsC1$FDR<0.05],rownames(C1vsC3)[C1vsC3$FDR<0.05],rownames(C2vsC3)[C2vsC3$FDR<0.05]) #consider only the statistically significant one
genes_of_interest=genes_of_interest[!duplicated(genes_of_interest)]
filtered_data <- expression_long %>% filter(gene %in% genes_of_interest)
filtered_data$expression=as.numeric(filtered_data$expression) #ensure numeric class

# Create the Violin Plot
custom_colors <- c("1" = "#FF6247" , "2" = "#3CB371", "3" = "#1E90FF")
ggplot(filtered_data , aes(x = Clusters, y = expression, fill = Clusters)) +
  geom_violin(trim = FALSE, alpha = 0.6) +  # Violin plot for distribution
  geom_boxplot(width = 0.1, outlier.shape = NA) +  # Option1:  Boxplot inside
 # geom_jitter(width = 0.2, size = 1, alpha = 0.5) +  # Option2: Individual points inside
  facet_wrap(~gene, scales = "free_y") +  # Separate plots per gene
  theme_minimal() +
  labs(title = "Violin Plot of Gene Expression by Cluster",
       x = "Cluster", y = "Expression Level") +
  scale_fill_manual(values = custom_colors)

