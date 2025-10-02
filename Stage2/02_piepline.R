#libraries
library(DESeq2)
library(pheatmap)

#Read and ignore lignes of commant
fc <- read.delim("gene_counts.txt", comment.char = "#", check.names = FALSE)
head(fc)

#count file
c_e_count <- read.delim('gene_counts.txt', header = T)
c_e_meta <- read.delim('metadata.tsv', header = T, stringsAsFactors = T)

#preview
head(c_e_count)
head(c_e_meta)

#keep to the important columns
raw_counts <- fc[, c("Geneid", grep("^SRR", colnames(fc), value = TRUE))]
head(raw_counts)

#add rownames
rownames(raw_counts) <- c_e_count$Geneid
head(raw_counts)
raw_counts <- raw_counts[ , -1]  # on enlève la colonne Geneid

#create desqdataset
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = c_e_meta,
                              design = ~State)

#preview
dds
dds$SRA.Accession.Number
dds$Age
dds$State

#show the ref
dds$State <- relevel(dds$State, ref = "Healthy_control")

design(dds) <- ~ Age + State

#Perform the differential expression analysis
dds <- DESeq(dds)

#with 0.05 we only had 3 genes, so in order to have more genes we choose 0.1
final_res <- results(dds, contrast = c("State", "Bipolar_disorder", "Healthy_control"), alpha = 0.10)

#look at your result
head(final_res)
summary(final_res)
sum(final_res$padj < 0.1, na.rm = TRUE)


#we have a truncated data, let's see the distro of p-values

plot(density(x = na.omit(final_res$pvalue)))

#ok let's look at our differentially expressed genes

res_tbl <- as.data.frame(final_res)

# Volcano plot de base

plot(x = final_res$log2FoldChange, 
     y = -log10(final_res$pvalue), 
     cex = 0.25,
     pch = 19, 
     col = 'grey',
     ylab = '-log10(p-value)',
     xlab = 'Log2 Fold Change')

abline(v = c(-1, 1), lty = 2)  # seuil LFC
abline(h = -log10(0.1), lty = 2)  # seuil p-value

# color points of volcano plot
upregulated <- subset(final_res, padj < 0.1 & log2FoldChange > 1)
downregulated <- subset(final_res, padj < 0.1 & log2FoldChange < -1)

points(upregulated$log2FoldChange,
       -log10(upregulated$pvalue),
       cex = 0.5, col = 'red', pch = 19)

points(downregulated$log2FoldChange,
       -log10(downregulated$pvalue),
       cex = 0.5, col = 'blue', pch = 19)

mtext('Volcano (FDR 0.10 & |LFC| > 1)')


#matrice transformed
vsd <- vst(dds, blind = FALSE)            # tient compte du design
mat <- assay(vsd)

#choose genes to show 2) choisir des gènes à afficher

sel <- with(as.data.frame(final_res),
            which(!is.na(padj) & padj < 0.10 & abs(log2FoldChange) > 1)) # DEGs reasonable

#if there is small amoutn of genes, take only the most variable ones
if (length(sel) < 30) {
  v <- rowVars(mat)
  sel <- head(order(v, decreasing = TRUE), 300)
}

mat_sub <- mat[sel, ]

#standarisation in order to compare the profiles 
mat_sub <- t(scale(t(mat_sub)))           # z-score per gene 

#columns annotation 
ann_col <- as.data.frame(colData(dds)[, c("State","Age")])
# colors
ann_colors <- list(
  State = c(Bipolar_disorder = "tomato", Healthy_control = "steelblue")
)

#heatmap
pheatmap(mat_sub,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         fontsize_col = 10,
         scale = "none")   


gene_list <- rownames(subset(final_res, padj < 0.1 & abs(log2FoldChange) > 1))
write.csv(gene_list, 'gene_list.csv')


#Functional Enrichment Analysis
# Visit https://bioinformatics.sdstate.edu/go/

