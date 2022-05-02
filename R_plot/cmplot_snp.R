library("CMplot")

# base_dir = "/home/galaxy/project/alleleSpecific_analysis/results/ASE_analysis/readCount/Fisher_test/by_ind/sig/cmplot/"
# base_dir = "/home/galaxy/project/alleleSpecific_analysis/results/ASm6A/cmplot/"
base_dir = "/home/galaxy/project/alleleSpecific_analysis/results/ASm6A_addZr/cmplot/"
setwd(base_dir)  # base_dir = "/home/galaxy/project/alleleSpecific_analysis/results/ASm6A/"

hela_cis = read.table("cm_format.txt", header = TRUE, sep = "\t")
CMplot(hela_cis, plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="pdf",memo="",dpi=600,
       file.output=TRUE,verbose=TRUE,width=9,height=6)

CMplot(hela_cis,plot.type="m",LOG10=FALSE,threshold=2,file="jpg",memo="",dpi=600,
       file.output=TRUE,verbose=TRUE,width=14,height=6,ylab="log2(SNP effect size)")
