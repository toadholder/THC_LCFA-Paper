#Install DESeq if you don't have it yet (from Bioconductor) 
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "EnhancedVolcano")) 

library(DESeq2)

#load in raw, unnormalized counts 
library(readxl)
GQ_rawcounts_withtopups <- read_excel("C:/Users/temea/OneDrive - Queen's University/Desktop/Paper/LCFA_THC Paper/GQ_rawcounts_withtopups.xlsx", 
                                      +     col_types = c("text", "text", "numeric", 
                                                          +         "numeric", "numeric", "numeric", 
                                                          +         "numeric", "numeric", "numeric", 
                                                          +         "numeric", "numeric", "numeric", 
                                                          +         "numeric", "numeric", "numeric", 
                                                          +         "numeric", "numeric", "numeric"))

#Format file to set row names as Ensembl gene ID (gene names preserved in original dataframe) 
library(tibble)
Counts<-GQ_rawcounts_withtopups[, -2]
rownames(Counts)<-Counts$gene_id
names<-Counts
Counts<-Counts[, -1]
rownames(Counts)<-names$gene_id
rm(names)

#create colData table; column names are sample number
#Samples 1-4 = "VEH", Samples 6-8, 21 = "VEH_O3", Samples 11-14 = "THC", Samples 16-18, 22 = "THC_O3"
condition<-factor(c(rep("VEH",4), rep("VEH_O3",4), rep("THC", 4), rep("THC_O3", 4)))
coldata<-data.frame(row.names=colnames(Counts), condition)

#Set up DESeq2 object 
dds <-DESeqDataSetFromMatrix(countData=Counts, colData=coldata, design=~condition)

#Filter to remove low count genes 
smallestGroupSize<-4
keep<-rowSums(counts(dds) >=10) >= smallestGroupSize 
dds<-dds[keep,]
#15 900 genes left after filtering 


#Run DESeq with design considering 4 groups (use Wald test first for pair-wise comparisons)
dds<-DESeq(dds, test = "Wald")

#Apply variance-stabilizing transformation to account for highly expressed genes being more variable, masking patterns in lowly expressed genes
#Use blind = False modifier to allow vst to account for groups (i.e., not blinded)
vsd<-vst(dds, blind = FALSE)

#PCA plot
plotPCA(vsd, intgroup = "condition")
#Note that this function is using ntop=500 top features by variance 
#Sample distances (Euclidean) and plot as heatmap
distances<-dist(t(assay(vsd)))
heatmap(as.matrix(distances))
#Interpreting: darker = samples more different in expression profile, lighter = more similar

#Use LRT instead of Wald test to run omnibus test ("Does expression differ between the treatment groups?")
dds_lrt<-DESeq(dds, test="LRT", reduced = ~1)
#Results 
res_lrt<-results(dds_lrt, independentFiltering=TRUE, alpha=0.05)
summary(res_lrt)
#Any significant genes with FDR of 5%?
sig<-subset(res_lrt, padj<0.05)
summary(sig)
#1792 genes with adjusted p-value of < 0.05, varying among any groups


#Create a matrix that highlights which comparison is driving the LRT identified (significant) gene's differential expression
sig_genes<-rownames(subset(res_lrt, padj<0.05))
groups<-c("VEH", "VEH_O3", "THC", "THC_O3")
pairs<-combn(groups, 2, simplify=FALSE) #This generates all combinations for the groups 
pairwise_LFC<-matrix(NA, nrow=length(sig_genes), ncol=length(pairs))
rownames(pairwise_LFC)<-sig_genes
colnames(pairwise_LFC)<-sapply(pairs, function(x) paste(x[2], "vs", x[1], sep="_"))
for (i in seq_along(pairs)) {
  comp<-pairs[[i]]
  res<-results(dds, contrast = c("condition", comp[2], comp[1]))
  pairwise_LFC[, i]<-res[sig_genes, "log2FoldChange"]
}

driver_comparison<-apply(pairwise_LFC, 1, function(x) {
  colnames(pairwise_LFC)[which.max(abs(x))]
})

driver_LFC<-apply(pairwise_LFC, 1, function(x) {
  x[which.max(abs(x))]
})

#Make summary table 
drivers<-data.frame(
  Gene = sig_genes, 
  Driving_Comparison=driver_comparison,
  LFC=driver_LFC
)
#Export this table 
write.csv(drivers, file = "driver_genes.csv", row.names=FALSE)

#Make heatmap to illustrate driver genes in context of treatment group comparisons
install.packages("pheatmap")
library(pheatmap)
#Order genes by driving comparison
annotation_row<-data.frame(Driving=driver_comparison)
rownames(annotation_row)<-rownames(pairwise_LFC)
gene_order<-order(annotation_row$Driving)
pairwise_LFC_ordered<-pairwise_LFC[gene_order, ]
annotation_row_ordered<-annotation_row[gene_order, , drop=FALSE]

pheatmap(
  pairwise_LFC_ordered,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  annotation_row = annotation_row_ordered,
  show_rownames = FALSE,
  main = "Pairwise LFCs for LRT-significant genes"
)
#Interpretation: Rows are the LRT significant genes grouped by pairwise comparison with largest log2 fold change
#Colours show magnitude and direction of effect; supposed to show which genes are driving differences among groups
#Perform hierarchical clustering of genes and label the clusters
hc<-hclust(dist(pairwise_LFC_ordered), method="ward.D2")
#Cut tree into clusters, can adjust k but use k=4
k<-4
clusters<-cutree(hc, k=k)
annotation_row_ordered$Cluster<-factor(clusters)

#Summarize driving comparisons per cluster (which driving comparison is most frequent for each cluster?)
library(dplyr)
cluster_summary <- annotation_row_ordered %>%
  as.data.frame() %>%
  group_by(Cluster, Driving) %>%
  summarise(Count = n(), .groups = "drop")

#Make heatmap 
pheatmap(
  pairwise_LFC_ordered,
  scale = "row",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  annotation_row = annotation_row_ordered,
  show_rownames = FALSE,
  main = "Pairwise LFCs for LRT-significant genes with Cluster Labels"
)

#For biological significance, generate table containing which genes are clustered/clustering patterns
#Data comes from the hierarchical clustering of LFCs for sig. genes identified by DEA (LRT method)
cluster_table<-data.frame(
  Gene=rownames(pairwise_LFC_ordered),
  Cluster=clusters
)

#Split into lists of genes per cluster
genes_by_cluster<-split(cluster_table$Gene, cluster_table$Cluster)

#Add in which treatment comparison is driving each gene 
cluster_table$Driving<-annotation_row_ordered$Driving

#Export as Excel file for reporting and downstream enrichment analysis reporting 
write.csv(cluster_table, "4Way_DEA_geneclusters.csv", row.names=FALSE)