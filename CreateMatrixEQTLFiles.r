## This script uses the output of the QC.r scripts and formats the methylation, genotype and covariate data for use in MatrixEQTL.
## Please consult the webpage for MatrixEQTL for further details on how to format data for this package: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
## Filepaths and sample IDs have been removed for security reasons, therefore it serves as a guide to the analysis only 

setwd("")
load("QCed_Data.rda")

crosshyb<-read.csv("/data/Minerva/FromEilis/450KAnno/CrossHybridisingProbesPriceORWeksberg.csv", row.names = 1)
probes<-read.csv("/data/Minerva/FromEilis/450KAnno/SNPsinProbesAnno.csv", row.names = 1)

## remove cross hybridising probes
remove<-match(crosshyb[,1], rownames(betas))
remove<-remove[which(is.na(remove) != TRUE)]
betas<-betas[-remove,]
## remove SNP probes as not a twin analysis
probes<-probes[row.names(betas),]
betas<-betas[which(probes$Weksburg_CommonSNP_Af_within10bpSBE == "" | probes$Illumina_CommonSNP_Af_within10bpSBE == ""),]
betas<-betas[-grep("rs", rownames(betas)),]


## Exclude sex chromosomes
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
probeAnnot<-probeAnnot[rownames(betas),]
probeAnnot<-probeAnnot[which(probeAnnot$chr != "chrX"),]
probeAnnot<-probeAnnot[which(probeAnnot$chr != "chrY"),]
betas<-betas[rownames(probeAnnot),]

## load geno data to match
geno<-read.table("_chr22.raw", header = T, stringsAsFactors = FALSE)
geno.map<-read.table("_chr22.bim", stringsAsFactors = FALSE)
colnames(geno.map)<-c("Chr", "id", "cm", "bp", "A1", "A2")

geno<-geno[match(intersect(geno$IID, sampleSheet$Sample_Name), geno$IID),]
betas<-betas[, match(geno$IID, sampleSheet$Sample_Name)]
sampleSheet<-sampleSheet[match(geno$IID, sampleSheet$Sample_Name),]

write.table(betas, "MatrixEQTL/Methylation.txt", sep = "\t", quote = FALSE)

pca<-read.table("*.PCA.eigenvec", stringsAsFactors = FALSE)
pca<-pca[match(geno$IID, pca[,2]),]

cov<-t(cbind(as.numeric(sampleSheet$gender.x.chr)-1, pca[,3], pca[,4], pca[,5], pca[,6], pca[,7]))
colnames(cov)<-sampleSheet$Sample_Name
rownames(cov)<-c("Sex", "PC1", "PC2", "PC3", "PC4", "PC5")
write.table(cov, "MatrixEQTL/Covariates_Sex_5PCs.txt", sep = "\t", quote = FALSE)

geno.ids<-geno$IID

geno<-geno[,-c(1:6)]
geno<-t(geno)
colnames(geno)<-sampleSheet$Sample_Name

## filter to exclude genetic variants with < 5 counts per genotype group
a<-NULL
for(i in 1:nrow(geno)){
	if(min(table(geno[i,])) < 5){
		a<-append(a, i)
		}
}
geno<-geno[-a,]
geno.map<-geno.map[-a,]
write.table(geno, "MatrixEQTL/Genotypes_MinGenoCount5_Chr22.txt", sep = "\t", quote = FALSE)
write.table(cbind(geno.map$id, geno.map$Chr, geno.map$bp), "MatrixEQTL/Genotypes_MapInfo_Chr22.txt", sep = "\t", quote = FALSE)

for(chr in 1:21){
	geno<-read.table(paste("_chr", chr, ".raw", sep = ""), header = T, stringsAsFactors = FALSE)
	geno.map<-read.table(paste("_chr", chr, ".bim", sep = ""), stringsAsFactors = FALSE)
	colnames(geno.map)<-c("Chr", "id", "cm", "bp", "A1", "A2")

	geno<-geno[match(geno.ids, geno$IID),]
	
	geno<-geno[,-c(1:6)]
	geno<-t(geno)
	colnames(geno)<-sampleSheet$Sample_Name


	a<-NULL
	for(i in 1:nrow(geno)){
		if(min(table(geno[i,])) < 5){
			a<-append(a, i)
			}
	}
	geno<-geno[-a,]
	geno.map<-geno.map[-a,]
	write.table(geno, paste("MatrixEQTL/Genotypes_MinGenoCount5_Chr", chr, ".txt", sep = ""), sep = "\t", quote = FALSE)
	write.table(cbind(geno.map$id, geno.map$Chr, geno.map$bp), paste("MatrixEQTL/Genotypes_MapInfo_Chr", chr, ".txt", sep = ""), sep = "\t", quote = FALSE)
	
}

