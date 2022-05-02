library(DESeq2)
setwd("/home/galaxy/project/circRNA_analysis/sailfish-circ/count_matrix")
data = read.table("reads_matrix.txt", header = TRUE, row.names = 2)
head(data)
coldata <- read.table("colData.txt", sep = "\t", header = TRUE, row.names = 1)
head(coldata)
all(rownames(coldata) %in% colnames(data))

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
# Note that we could have specified the coefficient or contrast we want to build a results table for, using either of the following equivalent commands:
res <- results(dds, contrast=c("condition","trt","untrt"))  # , name="condition_trt_vs_untrt"
# or to shrink log fold changes association with condition:
res_2 <- lfcShrink(dds, coef="condition_untrt_vs_trt", type="apeglm")
resOrdered <- res_2[order(res_2$pvalue),]
resSig <- subset(resOrdered, padj < 0.05)
write.csv(as.data.frame(resSig), file="result_DESeq2.csv", quote = FALSE)

#
plotMA(res_2, main="DESeq2", ylim=c(-2,2))



###############################################################
###############   featurecounts && deseq2   #################
###############################################################
library(DESeq2)
library(Rsubread)
library(ggplot2)
setwd("/home/galaxy/project/alleleSpecific_analysis/results/fetus_vs_adult/diff_expression/featureCounts/")
# sampleNames <- c("heart_1_L6", "heart_2_L7", "heart_3_L6", "CRR042291", "CRR055528")
# 第一行是命令信息，所以跳过
data <- read.table("muscle-GeneCounts.txt", header=TRUE, quote="\t", skip=1)
# 前六列分别是Geneid	Chr	Start	End	Strand	Length
# 我们要的是count数，所以从第七列开始
df <- data[,7:ncol(data)]
sampleNames <- colnames(df)
# names(data)[7:11] <- sampleNames
countData <- as.matrix(data[7:ncol(data)])
head(countData)
rownames(countData) <- data$Geneid
database <- data.frame(name=sampleNames, condition=c("fetus", "fetus", "adult", "adult","adult", "adult","adult"))
# database <- data.frame(name=sampleNames, condition=c("fetus","fetus","fetus","adult","adult","adult","adult","adult","adult","adult","adult","adult","adult","adult","adult"))
rownames(database) <- sampleNames
## 设置分组信息并构建dds对象
dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]
## 使用DESeq函数估计离散度，然后差异分析获得res对象
dds <- DESeq(dds)
res <- results(dds)
# res = res[order(res$pvalue),]
# diff_gene_deseq2 <- subset(res,padj<0.05 & abs(log2FoldChange)>1)
write.csv(res, "muscle_des_output.csv")
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata, "muscle_all_des_output.csv", row.names=FALSE)
#### 
### ID convertion
# "https://www.biotools.fr/human/refseq_symbol_converter"
################################## plot ##############################
############################# figure 1  #############################
plotMA(res, main="Heart DESeq2", ylim=c(-2,2))

############################# figure 2  #############################
# 这里的resdata也可以用res_des_output.csv这个结果重新导入哦。
# 现在就是用的前面做DESeq的时候的resdata。
resdata$change <- as.factor(
  ifelse(
    resdata$padj<0.01 & abs(resdata$log2FoldChange)>1,
    ifelse(resdata$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)
valcano <- ggplot(data=resdata, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=15) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) + 
  ggtitle("Heart DESeq2 Valcano") + 
  scale_color_manual(name="", values=c("red", "green", "black"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)

########
valcano
########

############################# figure 3  #############################
rld <- rlog(dds)
pcaData <- plotPCA(rld, intgroup=c("condition", "name"), returnData=T)
percentVar <- round(100*attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=name)) + 
  geom_point(size=3) + 
  ggtitle("DESeq2 PCA") + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))
pca


###################  bioMart ########################################
library('biomaRt')
library("curl")
mart <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
my_ensembl_gene_id<-row.names(diff_gene_deseq2)
#listAttributes(mart)
mms_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                    filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)

