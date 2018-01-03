## This script uses the output of the QC.r & methylationDerivedVariables.r scripts to perform the primary analysis of a genome-wide analysis to identify dna methylation sites associated with ASD status. 
## Filepaths and sample IDs have been removed for security reasons, therfore it serves as a guide to the analysis

uniqueAnno<-function(row){
if(row != ""){
	return(paste(unique(unlist(strsplit(row, "\\;"))), collapse = ";"))
	} else {
	return(row)
	}
	}
  
library(qqman)

setwd("")
load("QCed_Data.rda")

crosshyb<-read.csv("CrossHybridisingProbesPriceORWeksberg.csv", row.names = 1)
probes<-read.csv("SNPsinProbesAnno.csv", row.names = 1)

## remove cross hybridising probes
remove<-match(crosshyb[,1], rownames(betas))
remove<-remove[which(is.na(remove) != TRUE)]
betas<-betas[-remove,]
## remove SNP probes
probes<-probes[row.names(betas),]
betas<-betas[which(probes$Weksburg_CommonSNP_Af_within10bpSBE == "" | probes$Illumina_CommonSNP_Af_within10bpSBE == ""),]
betas<-betas[-grep("rs", rownames(betas)),]

## load cellular composition variables
counts<-read.csv("Analysis/HousemanEstimatedCellCounts.csv", stringsAsFactors = FALSE, row.names = 1)
counts<-counts[match(colnames(betas) ,rownames(counts)),]

### alternatively use smokingScores
scores<-read.csv("Analysis/QC/SmokingScores.csv", stringsAsFactors = FALSE, row.names = 1)
scores<-scores[match(colnames(betas), scores$Basename),]

pheno<-read.csv("", stringsAsFactors = FALSE)
### needs reformating
pheno.cont<-pheno[,c(22,18,grep("cont", colnames(pheno)))]
pheno.case<-pheno[,c(21,18,grep("case", colnames(pheno)))]
pheno.cont<-cbind(pheno.cont, "Con")
pheno.case<-cbind(pheno.case, "Case")
colnames(pheno.cont)<-gsub("_cont", "", colnames(pheno.cont))
colnames(pheno.case)<-gsub("_case", "", colnames(pheno.case))
colnames(pheno.cont)[1]<-"participiant_id"
colnames(pheno.case)[1]<-"participiant_id"
colnames(pheno.cont)[13]<-"CasCon"
colnames(pheno.case)[13]<-"CasCon"
pheno<-rbind(pheno.cont,pheno.case)
pheno$participiant_id<-gsub("MMXII_iPSYCH_", "", pheno$participiant_id)
pheno<-pheno[(match(sampleSheet$Sample_Name, pheno$participiant_id)),]

birthMonth<-unlist(strsplit(pheno$fdato, "/"))[seq(from = 2, by = 3, to = 3*nrow(pheno))]
birthMonth<-factor(birthMonth)
birthYear<-unlist(strsplit(pheno$fdato, "/"))[seq(from = 3, by = 3, to = 3*nrow(pheno))]
birthYear<-factor(birthYear)

res<-matrix(data = NA, nrow = nrow(betas), ncol = 3)
colnames(res)<-c("P-value", "MeanDiff", "SE")
for(i in 1:nrow(betas)){

	model<-lm(betas[i,] ~ as.factor(sampleSheet$CaCo) + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + birthMonth + birthYear + as.factor(pheno$geo5) + pheno$gest_week)
	res[i,]<-summary(model)$coefficients["as.factor(sampleSheet$CaCo)Ctrl", c(4,1,2)]
}

### as effects are increases in controls, multiple by -1
res[,2]<-res[,2]*-1

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
probeAnnot<-probeAnnot[rownames(betas),]
res<-cbind(res, probeAnnot)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other)
probeAnnot<-probeAnnot[rownames(betas),]
res<-cbind(res, probeAnnot)

res<-res[order(res[,1]),]

write.csv(res, "Analysis/CaseControl_Sex_Batch_HousemanCellComp_SmokingScore_BY_BM_GA_Urbanicity.csv")

pdf("Analysis/TopProbes_AUT_EWAS_HousemanCellComp_GA_Urbanicity.pdf", height = 6, width = 12)
par(mfrow = c(2,5))
par(mar = c(4.5,4.5,5,0.2))
for(i in 1:10){
	model<-lm(betas[rownames(res)[i],] ~ as.factor(sampleSheet$CaCo) + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + birthMonth + birthYear + as.factor(pheno$geo5) + pheno$gest_week)
	data.points<-(residuals(model)+coef(model)[2]*c(0,1)[as.factor(sampleSheet$CaCo)] + coef(model)[1])*100
	boxplot(data.points ~ as.factor(sampleSheet$CaCo), ylab = "DNA methylation (%)", main = paste(rownames(res)[i], "\n", uniqueAnno(res$UCSC_RefGene_Name[i]), sep = ""), col = c(rgb(1,0,0,0.25), rgb(0,1,0,0.25)), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
	title(main = paste("p = ", signif(res[i,"P.value"],3), sep = ""), line = 0, adj = 1)
}
dev.off()


res.all<-res.all[which(res.all$chr != "Y")]
res.all$chr<-gsub("chr", "", res.all$chr)
res.all$chr[which(res.all$chr == "X")]<-23
res.all$chr<-as.numeric(res.all$chr)
par(mar = c(5,5,2,1))
manhattan(res.all, p = "P.value", chr = "chr", bp = "pos", suggestiveline = -log10(5e-5), genomewideline = -log10(1e-7))
