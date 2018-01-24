## Run colocalisation analysis of genetic signals for DNA methylation and ASD following method in Giambartolomei et al.
## Tests ASD against all DNA methylation sites found within 250 kb of a ASD GWAS loci as identified by the PGC-iPSYCH meta-analysis.

removeAlleles<-function(text){
	tmp<-unlist(strsplit(unlist(strsplit(gsub("X", "", text), "\\.")), "_"))
	if(length(tmp) == 3){
		return(tmp)
	} else {
		return(c(tmp[1:2], paste(tmp[3:length(tmp)], collapse = "_")))
	}
}

library(coloc)
library(data.table)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

setwd("")

probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)

gwas<-fread("") ## Load ASD GWAS results
output<-NULL

regions<-read.table("Analysis/AdditionalAutismRegions.txt", header = TRUE, row.names = 1)
regions<-regions[1:3,] ## only first 3 are genome-wide significant

for(i in 1:3){

	chr<-regions$CHR[i]
	start<-regions$BP[i]-250000
	stop<-regions$BP[i]+250000
	expression_file_name<-paste("MatrixEQTL/Methylation_ProbesForASDColocAnalysis_Chr", chr, "_", start, "_", stop, ".txt", sep = "")
	SNP_file_name = paste("MatrixEQTL/Genotypes_MinGenoCount5_Chr", chr, ".txt", sep = "")
	snps_location_file_name = paste("MatrixEQTL/Genotypes_MapInfo_Chr", chr, ".txt", sep = "")
	
  ## load mQTL results
	mQTL<-read.table(paste("MatrixEQTL/Output/AllmQTLsforASDColocAnalysis_mQTL_chr", chr,"_", start, "_", stop, ".txt", sep = ""), stringsAsFactors = FALSE, header = TRUE)
	mQTL<-cbind(unlist(lapply(strsplit(mQTL$SNP, "_"), head, n = 1)), mQTL)
	
  ## load variants frequencies
	freq<-read.table("01_Genotypes/FilteredGenotypes/GT_hwe_mind_geno_maf_freq.frq", header = TRUE, stringsAsFactors = FALSE)
	gwas.sub<-gwas[which(gwas$CHR == chr),]
	
	mQTL<-cbind(mQTL, paste(freq$SNP, apply(apply(freq[,c("A1", "A2")], 1, sort), 2, paste, collapse = "_"), sep = "_")[match(mQTL[,7], freq$SNP)])
	colnames(mQTL)[ncol(mQTL)]<-"SNP2"
	mQTL$SNP2<-as.character(mQTL$SNP2)
	
	dataset2<-list(beta= log(gwas.sub$OR), varbeta= (gwas.sub$SE^2), type = "cc", snp = gwas.sub$SNP, MAF = (gwas.sub$FRQ_A_18381*18381+gwas.sub$FRQ_U_27969*27969)/(18381+27969), N=46350, s = 0.397)

	for(each in unique(mQTL$gene)){
		mQTL.sub<-mQTL[which(mQTL$gene == each),]
		if(nrow(mQTL.sub) > 1){
			dataset1<-list(beta=mQTL.sub$beta,varbeta=(mQTL.sub$beta/mQTL.sub$t.stat)^2,type = "quant",snp = mQTL.sub[,1], MAF = freq$MAF[match(mQTL.sub[,1], freq$SNP)], N=1157)
			my.res<-coloc.abf(dataset1, dataset2)
			output<-rbind(output, c("trait2"=each, my.res$summary, i))
			}
		}
	}
}


sum<-as.numeric(output[,6])+as.numeric(output[,7])
ratio<-as.numeric(output[,7])/as.numeric(output[,6])
output<-cbind(output, sum, ratio)

## Annotate results
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
probeAnnot<-probeAnnot[match(output[,1], rownames(probeAnnot)),c("chr", "pos")]
output<-cbind(output, probeAnnot)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other)
probeAnnot<-probeAnnot[match(output[,1], rownames(probeAnnot)),5:14]
output<-cbind(output, probeAnnot)

write.csv(output, "Coloc/ColocAnalysis_ASD_GenomeWideSigRegionsOnly.csv")

