#annotation
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(plyr)
library(pathview)

#setwd("/Data_Resource/m6A-seq/com_fetus_vs_adult/diff_expression/DESeq2/heart")
# setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/fetal_vs_adult/specific/filter/overlap_gene/by_tissue/")
# setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/by_tissue/overlap_gene/")
# setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/by_tissue/overlap_GWAS/overlap_gene/")
# setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/fetal_vs_adult/specific/filter/overlap_gene/by_tissue/merge/")
# setwd("/Charles/project/ASm6A/ASm6A/by_tissue/overlap_gwas/target_genes/")
# setwd("/Charles/project/ASm6A/ASm6A/by_tissue/target_genes/")
setwd("/Charles/project/ASm6A/ASm6A/by_tissue/overlap_gene/")
# GO annotation
go_anno <- function(gene_symbol, data_type){
  # genes <- gene_entrezid$SYMBOL
  for(an_type in c("BP", "MF", "CC")){
    print(paste("GO annotation_", an_type, sep = ""))
    ggo <- groupGO(gene = gene_symbol,
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = an_type,
                   level = 3)
    file_prefix <- paste("GO_ANNO_", an_type, "_", data_type, sep = "")
    write.table(ggo, paste(file_prefix, ".txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
    pdf(file = paste(file_prefix, ".pdf", sep = ""), width = 8, height = 8)
    print(barplot(ggo, drop=TRUE, showCategory=12))
    dev.off()
  }
}

# go_anno(pos_data, "positive")
# go_anno(neg_data, "negative")

# GO enrichment
go_enri <- function(gene_symbol, data_type){
  # genes <- gene_entrezid$SYMBOL
  for(an_type in c("BP", "MF", "CC")){
    print(paste("GO Enrichment_", an_type, sep = ""))
    ego <- enrichGO(gene = gene_symbol,
                    OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL",
                    ont = an_type,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
    file_prefix <- paste("GO_ENRICH_", an_type, "_", data_type, sep = "")
    write.table(ego, paste(file_prefix, ".txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
    pdf(file = paste(file_prefix, ".pdf", sep = ""), width = 8, height = 8)
    print(barplot(ego)) # , showCategory=12
    dev.off()
  }
}
#
# select significant
#sub_ego <- subset(as.data.frame(ego), p.adjust < 0.05 & qvalue < 0.05)
#sub_ego <- ldply(sub_ego, data.frame)



# go_enri(pos_data, "positive")
# go_enri(neg_data, "negative")

############################################################################
# KEGG
###########################################################################
# KEGG annotation
# KEGG Enrichment
# transform SYMBOL to ENTREZID
# pos_gene_id <- bitr(pos_data$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
# neg_gene_id <- bitr(neg_data$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")$ENTREZID
#
kegg_rich <- function(gene_entrezid, data_type){
  print("KEGG Enrichment")
  kk <- enrichKEGG(gene = gene_entrezid,
                   keyType = 'ncbi-geneid',
                   organism = "hsa",
                   pvalueCutoff = 0.07)
  # print(head(kk))
  file_prefix <- paste("KEGG_ENRICH_", data_type, sep = "")
  write.table(kk, paste(file_prefix, ".txt", sep = ""), row.names = FALSE, quote = FALSE, sep = "\t")
  pdf(file = paste(file_prefix, ".pdf", sep = ""), width = 8, height = 8)
  print(barplot(kk, showCategory = 14))
  dev.off()
  # plot each kegg id pathway
  #kegg_id <- summary(kk)$ID
  #for(i_id in kegg_id){
  #pdf(file = paste("KEGG_pathway_", i_id, ".pdf", sep = ""))
  #print(browseKEGG(kk, i_id))
  #dev.off()
  #}
}

# kegg_rich(pos_gene_id, "positive")
# kegg_rich(neg_gene_id, "negative")


reactome_rich <- function(gene_entrezid, data_type){
  print("Reactome Enrichment")
  kk <- enrichPathway(gene = gene_entrezid,
                   readable=T,
                   pvalueCutoff = 0.05)
  # print(head(kk))
  file_prefix <- paste("Reactome_ENRICH_", data_type, sep = "")
  pdf(file = paste(file_prefix, ".pdf", sep = ""), width = 14, height = 8)
  print(barplot(kk, showCategory = 12))
  dev.off()
}


for (file in list.files(path="./",pattern="*.bed")){
  print(file)
  data <- read.table(file, sep = "\t", header = FALSE) # , header = TRUE
  # gene_symbol <- data$gene
  gene_symbol <- data$V8 # v4
  gene_symbol <- as.character(gene_symbol)
  # translate_results = bitr(gene_refseq, fromType = "REFSEQ", toType = c("SYMBOL", "ENTREZID"), OrgDb = "org.Hs.eg.db")
  translate_results = bitr(gene_symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  #gene_symbol = translate_results$SYMBOL
  gene_entrezid = translate_results$ENTREZID
  # print(head(gene_entrezid))
  cluster_name <- strsplit(file, ".bed")[[1]][1]
  print(cluster_name)
  #go_anno(gene_symbol, cluster_name)
  #go_enri(gene_symbol, cluster_name)
  kegg_rich(gene_entrezid, cluster_name)
  # reactome_rich(gene_entrezid, cluster_name)
}

# 检验是否gene在kegg中没有注释；
# https://guangchuangyu.github.io/2017/11/kegg-gene-annotation/
# bitr_kegg(geneID = gene_entrezid, fromType='ncbi-geneid', toType='Path', organism='hsa') -> xx
