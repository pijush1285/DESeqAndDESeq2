


library(DESeq2)

#This is a veru good website to do differential expression analysis of RNA seq data.
#http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

#####################################################################
# Human Samples
countData <- as.matrix(read.csv("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountHuman/Rcount_Human/gene_count_matrix.csv", row.names="gene_id"))
# Mouse Samples
countData <- as.matrix(read.csv("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountMouse/Count/gene_count_matrix.csv", row.names="gene_id"))
countData <- as.matrix(read.csv("/home/pijush/Desktop/NGC/NewTask4institute_21.06.21/SINP/ReadCount/ReadCountHuman/Rcount_Human/transcript_count_matrix.csv", row.names="transcript_id"))

#####################################################################
# Second run using only the protein coding genes
# Human Samples
countData <- as.matrix(read.csv("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountHuman_v1/RCountH/Homo_sapiens.GRCh38.104.geneonly.csv", row.names="gene_id"))

# Mouse Samples
countData <- as.matrix(read.csv("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountMouse_v1/RCountM/Mus_musculus.GRCm39.104.geneonly.csv", row.names="gene_id"))
#####################################################################

#We need to compare :
#  A1 with B1.
countData1<-countData[,-c(2,4)]

#We need to compare :
#  A2 with B2.
countData2<-countData[,-c(1,3)]


#create experiment labels (two conditions)
colData <- DataFrame(condition=factor(c("A","A","B","B")))

#  A1 with B1.
colData <- DataFrame(condition=factor(c("A","B")))


# create DESeq input matrix
dds <- DESeqDataSetFromMatrix(countData, colData,formula(~ condition))
dds <- DESeqDataSetFromMatrix(countData1, colData,formula(~ condition))
# run DEseq
dds <- DESeq(dds)

# get deferentially expressed genes 
res <- results(dds)
View(res)
# order by BH adjusted p-value
resOrdered <- res[order(res$padj),]
# top of ordered matrix
head(res)
#To write the output.
setwd("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountHuman_v1/DifferentialExpressionAnalysis")
write.csv(res, file="Human_Diff_exp_Deseq2.gene_count.csv")
write.csv(res, file="P42vsP47_Diff_exp_Deseq2.transcript_count.csv")


setwd("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountMouse_v1/DifferentialExpressionAnalysis")
write.csv(res, file="Mouse_Diff_exp_Deseq2.gene_count.csv")

#############################################################################
# Valcano plots
#################

library("magrittr")
library("EnhancedVolcano")


EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'A versus B',
                subtitle = "Differential expression"
               )

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'N061011 versus N61311',
                pCutoff = 10e-32,
                FCcutoff = 0.5,
                pointSize = 3.0,
                labSize = 6.0)


############################################################################

res1 <- data.frame(res)
res1 <- data.matrix(res)
sigGenes <- res1[res1$padj <= 0.05 & abs(res1$log2FoldChange) > 2,]
sigGenes <- res1[abs(res1$log2FoldChange) > 2,]
sigGenes <- rownames(res1[res1$padj <= 0.05 & abs(res1$log2FoldChange) > 2,])

# Keep only the significantly differentiated genes where the fold-change was at least 3
sigGenes <- rownames(res[res$padj <= 0.05 & abs(res$log2FoldChange) > 1,])
sigGenes <- rownames(res[res$padj <= 0.05,])
sigGenes <- res[res$padj <= 0.1,]

res1 <- results(dds,alpha=0.05, lfcThreshold=1, altHypothesis="greaterAbs")
head(res1)
res1 <- res[abs(res$log2FoldChange) > 1, ]

res1 <- results(dds,alpha=0.05, lfcThreshold=1, altHypothesis="greaterAbs",independentFiltering = TRUE)

write.csv(res1, file="P42vsP47_Diff_exp_Deseq2.gene_count.csv")

res05 <- results(dds, alpha=0.05)
write.csv(res05, file="res05_Diff_exp_Deseq2.gene_count.csv")

#############################################################################

library("pheatmap")

select <- order(rowMeans(counts(dds, normalized=TRUE)),decreasing=TRUE)[1:20]
df<- as.data.frame(colData(dds)[ ,c("condition")])
ntd <- normTransform(dds)
pheatmap(assay(ntd)[select,], cluster_rows = FALSE, 
         show_rownames = FALSE, cluster_cols = FALSE, annotation_col = df)

pheatmap(assay(ntd)[select,])

###########################################################################

Diff_exp_Deseq2_gene_count_Excel <- read_excel("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountHuman/Differential Exp/Diff_exp_Deseq2.gene_count_Excel.xlsx", 
                                                 sheet = "fc 1.5, adjp 0.05")
View(Diff_exp_Deseq2_gene_count_Excel)
Diff_exp_Deseq2_gene_count_Excel$Gene_ID[1:10]

ntd <- normTransform(dds)
x <- assay(ntd)[Diff_exp_Deseq2_gene_count_Excel$Gene_ID,]
pheatmap(x)





library(readxl)
Diff_exp_Deseq2_gene_count_Excel <- read_excel("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountMouse/Differential Exp/Diff_exp_Deseq2.gene_count_Excel.xlsx", 
                                                 sheet = "FC 1.5, Adjp less 0.05")
View(Diff_exp_Deseq2_gene_count_Excel)

ntd <- normTransform(dds)
x <- assay(ntd)[Diff_exp_Deseq2_gene_count_Excel$Gene_ID,]
pheatmap(x)

#######################################################################
# Heat map 2nd times

library("pheatmap")
library(readxl)
Selected_Diffexpress_Genes_Adjp_less0_05andFC1_5 <- read_excel("ReadCountHuman_v1/Heatmap/Selected_Diffexpress_Genes_Adjp_less0.05andFC1.5.xlsx")
View(Selected_Diffexpress_Genes_Adjp_less0_05andFC1_5)

Selected_Diffexpress_Genes_Adjp_less0_05andFC1_5$Gene_ID[1:10]

ntd <- normTransform(dds)
x <- assay(ntd)[Selected_Diffexpress_Genes_Adjp_less0_05andFC1_5$Gene_ID,]
pheatmap(x)


###########
#Split row name only keep the gene name.

library(dplyr)
library(stringr)

str<-rownames(x)[1]
#Split the string only 
unlist(strsplit(str, "\\|"))
#Split string and grab the second one
unlist(strsplit(str, "\\|"))[2]


GeneList=NULL
str<-rownames(x)
for (i in 1:length(str)) {
  
  #If found no gene name then write only the ensg name.
  if (is.na(unlist(strsplit(str[i], "\\|"))[2])){
    print("First one Grab.") 
    print(str[i])
    GeneList[i]<-str[i] 
  } else
  {
    print("Second one grab.")
    print(unlist(strsplit(str[i], "\\|"))[2])
    GeneList[i]<-unlist(strsplit(str[i], "\\|"))[2]
  }
}

#Chaanging the row name.
rownames(x) <- GeneList
#Now Heatmap 
pheatmap(x)


########################################################################
#StatQuest: PCA in R
countData <- as.matrix(read.csv("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountHuman/Rcount_Human/gene_count_matrix.csv", row.names="gene_id"))
countData <- as.matrix(read.csv("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountMouse/Count/gene_count_matrix.csv", row.names="gene_id"))
head(countData)
countData[countData==0]<-1
pca <- prcomp(t(countData), scale. = TRUE) 

pca <- prcomp(t(countData)) 
plot(pca$x[,1], pca$x[,2])

setwd("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountHuman/PCA")
write.csv(pca$x, "Human_PCA_matrix.csv")

setwd("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountMouse/PCA")
write.csv(pca$x, "Mouse_PCA_matrix.csv")


library(ggfortify)
autoplot(pca)
autoplot(prcomp(t(countData)), label=TRUE, label.size=3, colour = c('red', 'green', 'blue', 'black'))
autoplot(prcomp(t(countData)), label=TRUE, label.size=3, colour = c('red', 'red', 'blue', 'blue'))





##################################################################################################
library("readr")
FPKM <- as.matrix(read_csv("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountHuman/FPKM_Human/FPKM Table.csv", row.names="gene_id"))
FPKM <- as.matrix(read_csv("C:/Users/PIJHUSH/Desktop/ReadCount SINP/ReadCountHuman/FPKM_Human/FPKM Table.csv"))
rownames(FPKM) <- FPKM[,1]
FPKM <- FPKM[,-1]
FPKM1<- data.frame(FPKM)
pca <- prcomp(t(FPKM1)) 




###############################################################
#Test with some MA plot 

res.noshr <- results(dds, name="dex_trt_vs_untrt")
plotMA(res.noshr, ylim = c(-5, 5))

plotMA(res, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})


hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")


