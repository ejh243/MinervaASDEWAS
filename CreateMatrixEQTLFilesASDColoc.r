## For colocalisation analysis need mQTL P values for all SNPs not just significant ones
## This script extracts DNA methylation sites within GWAS loci and creates files for MatrixEQTL.

setwd("")
meth<-read.table("MatrixEQTL/Methylation.txt", header = TRUE, row.names = 1)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
probeAnnot<-probeAnnot[rownames(meth),]


### only first 3 are genome-wide significant
regions<-read.table("Analysis/AdditionalAutismRegions.txt", header = TRUE, row.names = 1)


for(i in 1:3){

chr<-regions$CHR[i]
start<-regions$BP[i]-250000
stop<-regions$BP[i]+250000
probes<-as.character(rownames(probeAnnot)[which(probeAnnot$chr == paste("chr", chr, sep = "") & 
  probeAnnot$pos <= stop & probeAnnot$pos >= start)])
write.table(meth[intersect(probes, rownames(meth)),], 
  paste("MatrixEQTL/Methylation_ProbesForASDColocAnalysis_Chr", chr, "_", start, "_", stop, ".txt", sep = ""), 
  sep = "\t", quote = FALSE)

}
