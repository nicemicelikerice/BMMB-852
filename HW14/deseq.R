#Import raw read files and combine them
SEEM_WT1 <- read.delim("03-readcounts/SEEM_WT1_batch.txt", header=FALSE)
SEEM_WT2 <- read.delim("03-readcounts/SEEM_WT2_batch.txt", header=FALSE)
SEEM_WT3 <- read.delim("03-readcounts/SEEM_WT3_batch.txt", header=FALSE)

SEDM_WT1 <- read.delim("03-readcounts/SEDM_WT1_batch.txt", header=FALSE)
SEDM_WT2 <- read.delim("03-readcounts/SEDM_WT2_batch.txt", header=FALSE)
SEDM_WT3 <- read.delim("03-readcounts/SEDM_WT3_batch.txt", header=FALSE)

merge <-cbind(SEEM_WT1,SEEM_WT2[2],SEEM_WT3[2],SEDM_WT1[2],SEDM_WT2[2],SEDM_WT3[2])

row.names(merge) <- merge[,1]
merge <- merge[,-1]
colnames(merge) <- c("SEEM_WT1","SEEM_WT2","SEEM_WT3","SEDM_WT1","SEDM_WT2","SEDM_WT3")

merge <- head (merge,-5)

mydesign <- read.csv("design.csv", row.names=1)

write.csv(as.data.frame(merge),file="results/readcounts.csv")


#DESeq2
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = merge,colData = mydesign,design = ~ Time)
design(dds)<-formula(~Time)
dds$group <- factor(paste0(dds$Time))
design(dds)<-~group
dds <- DESeq(dds)
resultsNames(dds)


#Merge with normalized count data
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

#transformed read count
rld<-rlog(dds,blind=FALSE)
head(assay(rld),3)
write.csv(assay(rld),file="results/normalized data.csv")

vsd <- vst(dds, blind=FALSE)

#DEG_SEEM vs. SEDM
EvD <- results (dds, contrast=c("group","Five","Six"))
write.csv(as.data.frame(EvD),file="results/EvD_deseq2.csv")

#PCA plot using ggplot2
library(ggplot2)
data <- plotPCA(rld,intgroup=c("Time"),returnData=TRUE)
percentVar <- round(100*attr(data,"percentVar"))
pca_plot <- ggplot(data,aes(PC1,PC2,color=Time))+geom_point(size=10)+xlab(paste0("PC1:",percentVar[1],"%variance"))+ylab(paste0("PC1:",percentVar[2],"%variance"))
ggsave("plots/PCA_plot.png", plot = pca_plot, width = 8, height = 6, dpi = 300)

#MA plot
png("plots/MAplot.png",units="in",width=10,height=10,res=350)
plotMA(EvD,alpha = 0.1, main = "EvD", ylim=c(-10,10))
dev.off()

#Heatmap
library("pheatmap")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 30)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd))
anno$group <- NULL
pheatmap(mat, annotation_col = anno, main = "EvD heatmap", filename = "plots/heatmap.png")
