# setwd("/home/galaxy/project/alleleSpecific_analysis/shell/ASm6A_analysis/")
library(motifbreakR)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(MotifDb)
# genomes <- available.genomes()
# available.SNPs()

#setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/common_ASm6A/")
#snp <- "common_ASm6A_forR.bed"
setwd("/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/motifbreakR/")
snp <- "total_forR.bed"
snps.mb <- snps.from.file(file = snp,
                          search.genome = BSgenome.Hsapiens.UCSC.hg19,
                          format = "bed")
#
pwm_dir = "/home/galaxy/project/alleleSpecific_analysis/data/AS/oRNAment/PWMs/rename/"
file_list = list.files(path = pwm_dir, pattern="*.PWM", full.names = TRUE)
temp = list.files(path = pwm_dir, pattern="*.PWM")
#### construction ####
motifMetadata <- DataFrame(providerName = temp,
                           providerId = temp,
                           dataSource = "oRNAment",
                           geneSymbol = "RBP",
                           geneId =NA, geneIdType = NA, proteinId = NA, proteinIdType = NA, 
                           organism = "Hsapiens",
                           sequenceCount = NA, #  
                           bindingSequence = NA, bindingDomain = NA, tfFamily = NA, experimentType = "high-throughput methods", 
                           pubmedID = NA)
read_file <- function(in_file){
  df <- t(read.table(in_file, header = T, row.names = 1))
  ## must be T, can't be U
  rownames(df) <- c("A","C","G","T")
  return(df)
}

motif_ppm = lapply(file_list, read_file)
names(motif_ppm) <- temp
# names(motif_ppm) <- paste(motifMetadata$organism,motifMetadata$dataSource,motifMetadata$providerName,sep = "-")
MotifList_RBP <- new("MotifList", elementMetadata = motifMetadata, manuallyCuratedGeneMotifAssociationTable = data.frame(), listData = motif_ppm)
####################################################################
batchResult <- motifbreakR(snpList = snps.mb, filterp = TRUE,
                           pwmList = MotifList_RBP,
                           threshold = 1e-4,
                           method = "ic",
                           bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                           BPPARAM = BiocParallel::MulticoreParam(workers = 25))

# plotMB(results = results, rsid = "rs1006140", effect = "strong")
res = "motifbreakR_results.txt"
# res = "motifbreakR_commonASm6A_results.txt"
write.table(batchResult, file =res, sep = "\t", quote = FALSE)
