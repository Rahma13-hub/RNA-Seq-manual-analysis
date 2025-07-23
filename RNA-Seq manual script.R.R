#needed packages
library(matrixTests)
library(genefilter)
library(ggplot2)
library(dplyr)
library(pheatmap)
#loading the gene expression data for tumor and normal data
tumor = as.matrix(read.csv("lusc-rsem-fpkm-tcga-t_paired.csv", row.names = 1))
normal = as.matrix(read.csv("lusc-rsem-fpkm-tcga_paired.csv", row.names = 1))
dim(tumor)
dim(normal)
#bind the two cases together
data = cbind(tumor, normal)
#explore the dimension of the whole data
dim(data)
#explore if there is any missing expression values
sum(is.null(data))
is.na(data)
#explore the data distribution using histogram plot
hist(data, col = "red")
#scaling the data using log2 transformation
hist(log2(data+1), col = "red")
##filter low count genes which have row mean lower than 1
data = data[rowMeans(data) > 1,]
#calculate the logged mean for each group
tum.mean = apply((log2(data+1))[,1:50], 1,  mean)
norm.mean = apply((log2(data+1))[,51:dim(data)[2]], 1, mean)
#calculate the fold change by taking the difference between the two means
fold = tum.mean-norm.mean
#visualize the distribution of the fold change
hist(fold, col = "red")
#doing the differential expression statistical tests
phenotype = as.data.frame(factor(rep(c("tumor", "normal"), c(50, 50))))
colnames(phenotype)= "group"
#using non parametric test as wilcoxon test if the data isn't normally distributed
w = row_wilcoxon_twosample(data[,1:50], data[,51:dim(data)[2]])
#correct the wilcoxon test p value using the FDR method
p.adjusted = p.adjust(w$pvalue, method = "fdr") 
#save the result in data frame contain the fold change and the p adjusted value
result = as.data.frame(cbind(fold, p.adjusted))
#chose the statistical significant (DEGs) based on the p adjusted value less than 0.05
#and the biological significant fold change more than 2
result.deg = result[result$p.adjusted < 0.05 & abs(result$fold)>2,]
#export the Degs for further analysis
write.csv(as.matrix(result.deg), file = "result.deg.csv", quote = F, row.names = T)
#visualization using volcano plot
result$significant = ifelse(result$p.adjusted < 0.05 & abs(result$fold) > 2,
                            "significant", "not significant")
result$logP = -log10(result$p.adjusted)
ggplot(result, aes(x = fold, y = logP))+
  geom_point(aes(colour = significant))+
  scale_color_manual(values = c("significant" = "red", "not significant" = "gray"))+
  geom_vline(xintercept = c(-2, 2), 
             linetype = "dashed", color = "blue")+
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "blue")+
  labs(title = "volcano plot", x = "log2 fold change", y = "-log10 adjusted p-value")+
  theme_minimal()
top.significant = result.deg[order(result.deg$p.adjusted), ][1:50, ]
write.csv(top.significant, "top.50.genes.csv", quote = FALSE, row.names = TRUE)
df <- data.frame(
  term = c("Glycoprotein", "N-linked glycosylation", "Signal", 
           "Cell-matrix adhesion", "Extracellular space", 
           "Plasma membrane", "Heparin binding", "Secreted", 
           "Anchoring junction", "Tetranectin domain", "Cell junction", 
           "Blood coagulation", "Vessel maturation", 
           "Plasminogen activation", "Extracellular region"),
  pValue = c(1.2e-05, 1.4e-04, 0.00132, 0.00152, 0.00405, 0.00454, 
             0.0065, 0.01069, 0.0108, 0.01205, 0.01625, 0.01672, 
             0.01672, 0.01672, 0.0175)
)
df$logP <- -log10(df$pValue)
#visualize the top 15 enriched terms using ggplot
ggplot(df, aes(x = reorder(term, logP), y = logP)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 15 Enriched Terms", x = NULL, y = "-log10(P-value)") +
  theme_minimal()
top_genes <- rownames(top.significant)
top_data <- data[top_genes, ]
log_data <- log2(top_data + 1)
rownames(phenotype) <- colnames(log_data)
anno_colors <- list(
  group = c(tumor = "firebrick3", normal = "steelblue")
)
#visualize of gene expression using heatmap
pheatmap(log_data,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = phenotype,
         annotation_colors = anno_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 7,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Heatmap of Top 50 DEGs (Manual RNA-seq Analysis)")
