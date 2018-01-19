## This script uses the output of the QC.r & methylationDerivedVariables.r scripts and extends our ASD case control analysis to include an interaction with gender
## Filepaths and sample IDs have been removed for security reasons, therefore it serves as a guide to the analysis only 

setwd("")
load("QCed_Data.rda")

crosshyb<-read.csv("/data/Non_iPSYCH_Projects/Minerva/FromEilis/450KAnno/CrossHybridisingProbesPriceORWeksberg.csv", row.names = 1)
probes<-read.csv("/data/Non_iPSYCH_Projects/Minerva/FromEilis/450KAnno/SNPsinProbesAnno.csv", row.names = 1)

## remove cross hybridising probes
remove<-match(crosshyb[,1], rownames(betas))
remove<-remove[which(is.na(remove) != TRUE)]
betas<-betas[-remove,]
## remove SNP probes as not a twin analysis
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
### reformat to match DNA methylation data
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

res<-matrix(data = NA, nrow = nrow(betas), ncol = 15)
colnames(res)<-c("Aut:Females:P-value", "Aut:Females:MeanDiff", "Aut:Females:SE","Aut:Males:P-value", "Aut:Males:MeanDiff", "Aut:Males:SE","Aut:P-value", "Aut:MeanDiff", "Aut:SE", "Gender:P-value", "Gender:MeanDiff", "Gender:SE", "Int:P-value", "Int:MeanDiff", "Int:SE")
rownames(res)<-rownames(betas)

for(i in 1:nrow(betas)){
	model<-lm(betas[i,] ~ as.factor(sampleSheet$CaCo) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + birthMonth + birthYear + as.factor(pheno$geo5) + pheno$gest_week, subset = which(sampleSheet$gender.x.chr == "F"))
	res[i,1:3]<-summary(model)$coefficients["as.factor(sampleSheet$CaCo)Ctrl", c(4,1,2)]

	model<-lm(betas[i,] ~ as.factor(sampleSheet$CaCo)  + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + birthMonth + birthYear + as.factor(pheno$geo5) + pheno$gest_week, subset = which(sampleSheet$gender.x.chr == "M"))
	res[i,4:6]<-summary(model)$coefficients["as.factor(sampleSheet$CaCo)Ctrl", c(4,1,2)]
	
	model<-lm(betas[i,] ~ as.factor(sampleSheet$CaCo) + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$CaCo)*as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + birthMonth + birthYear + as.factor(pheno$geo5) + pheno$gest_week)
	res[i,13:15]<-summary(model)$coefficients["as.factor(sampleSheet$CaCo)Ctrl:as.factor(sampleSheet$gender.x.chr)M", c(4,1,2)]
	res[i,7:9]<-summary(model)$coefficients["as.factor(sampleSheet$CaCo)Ctrl", c(4,1,2)]
	res[i,10:12]<-summary(model)$coefficients["as.factor(sampleSheet$gender.x.chr)M", c(4,1,2)]
}

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
probeAnnot<-probeAnnot[rownames(betas),]
res<-cbind(res, probeAnnot)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other)
probeAnnot<-probeAnnot[rownames(betas),]
res<-cbind(res, probeAnnot)

write.csv(res, "Analysis/Aut_InteractionwithGender_wCovariates.csv")
