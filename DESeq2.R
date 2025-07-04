library(DESeq2)
setwd(dir = "D:/Workshop_RNA/")
counts <- read.table("counts.txt", header = TRUE, row.names = 1, comment.char = '#')
counts <- counts[, 6:ncol(counts)]
View(counts)
colnames(counts) <- c("T0_1","T0_2","T0_3","T150_1","T150_2","T150_3")

coldata <- data.frame(
  row.names = c("T0_1","T0_2","T0_3","T150_1","T150_2","T150_3"),
  condition = c("control","control","control","treated","treated","treated")
)
coldata

#Check colnames and rownames
all(rownames(coldata) %in% colnames(counts))

##Create DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
# Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)


##Perform PCA
# Perform PCA on log-transformed normalized counts
log_normalized_counts <- log1p(normalized_counts)  # Log-transform to stabilize variance
pca <- prcomp(t(log_normalized_counts))  # Transpose to make samples rows

# Inspect PCA results
summary(pca)

# Load required libraries
library(ggplot2)
library(ggrepel)  # For better label placement

# Extract principal components and add metadata
pca_data <- as.data.frame(pca$x)  # PCA results
pca_data$condition <- coldata$condition  # Add conditions
pca_data$sample <- rownames(pca_data)  # Add sample names as labels

# Plot PCA with labels
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 4) +  # Plot points
  geom_text_repel(size = 3, max.overlaps = 10) +  # Add labels with better spacing
  theme_minimal() +
  labs(title = "PCA of Normalized Counts with Labels", x = "PC1", y = "PC2")

# Variance explained by each principal component
variance <- pca$sdev^2 / sum(pca$sdev^2) * 100
barplot(variance, main = "Variance Explained by Principal Components", xlab = "Principal Components", ylab = "Percentage of Variance")

###Continue DESeq
dds <- DESeq(dds)
res <- results(dds)
summary(res)
#sig_res <- res[which(res$padj < 0.05), ]
sig_res <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 2), ]
summary(sig_res)

##Save output
write.csv(sig_res, "differential_expression_results.csv")

#Plot MA
plotMA(res, ylim=c(-2, 2))

# After running DESeq2
plotDispEsts(dds, main="Dispersion Estimates")

###Plot Volcano 
# Load required library
library(ggplot2)

# Define thresholds for significance
log2FC_cutoff <- 2  # Log2 fold change threshold
padj_cutoff <- 0.05 # Adjusted p-value threshold

# Add a column to classify significance
res$Significance <- "Not Significant"
res$Significance[res$padj < padj_cutoff & abs(res$log2FoldChange) > log2FC_cutoff] <- "Significant"

# Create the volcano plot
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +  # Add points
  scale_color_manual(values = c("gray", "red")) +  # Color mapping
  theme_minimal() +  # Clean theme
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value",
    color = "Significance"
  ) +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "blue") +  # FC threshold lines
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "blue")  # P-value threshold line


##Label Specific Genes
library(ggrepel)
ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value", color = "Significance") +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "blue") +
  geom_text_repel(data = subset(res, padj < padj_cutoff & abs(log2FoldChange) > log2FC_cutoff),
                  aes(label = rownames(subset(res, padj < padj_cutoff & abs(log2FoldChange) > log2FC_cutoff))),
                  size = 3)

##Heatmap visualization
library(pheatmap)
pheatmap(assay(rlog(dds))[rownames(sig_res), ])

# Transform counts for visualization
rld <- rlog(dds)  # Regularized log transformation
top_genes <- head(rownames(res[order(res$padj), ]), 20)  # Top 20 DE genes
pheatmap(assay(rld)[top_genes, ], cluster_rows=TRUE, cluster_cols=TRUE)


##Clusterprofiler
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

library(KEGGREST)
kegg_organisms <- keggList("organism")
head(kegg_organisms)  # Show the first few entries

organism <- "sce"
gene_list <- rownames(sig_res)
gene_list

# Perform KEGG enrichment analysis
kegg_enrich <- enrichKEGG(
  gene         = gene_list,
  organism     = organism,
  pvalueCutoff = 0.05
)
# Visualize results
barplot(kegg_enrich, showCategory = 10)
dotplot(kegg_enrich, showCategory = 10)
write.csv(as.data.frame(kegg_enrich), "kegg_enrichment_results.csv")

# Perform GO Enrichment Analysis
# Load the yeast OrgDb package
BiocManager::install("org.Sc.sgd.db")
library(org.Sc.sgd.db)
library(AnnotationDbi)
keytypes(org.Sc.sgd.db)
mapped_genes <- mapIds(org.Sc.sgd.db,keys = gene_list, column = "SGD",keytype = "ENSEMBL")
print(mapped_genes)
# Perform GO enrichment analysis
go_enrich <- enrichGO(
  gene         = mapped_genes,
  OrgDb        = org.Sc.sgd.db,
  keyType      = "SGD",         # Key type for yeast genes
  ont          = "BP",          # Biological Process
  pvalueCutoff = 0.05
)

# View results
head(go_enrich)

# Visualize results
barplot(go_enrich, showCategory = 10)
dotplot(go_enrich, showCategory = 10)

write.csv(as.data.frame(go_enrich), "go_enrichment_results.csv")
