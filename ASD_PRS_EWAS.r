## This script was used to perform an EWAS of ASD polygenic risk score. The polygenic risk scores were calculated in plink.
## Filepaths and sample IDs have been removed for security reasons, therfore it serves as a guide to the analysis

library(corrplot)
library("PerformanceAnalytics")

setwd("")

scores.all<-read.table("PRS/scores_iPSYCHwoMINERvA_AllSNPs/iPSYCHwoMINERvA_AllSNPs.S1.profile", header = TRUE)
for(i in 2:10){
	scores.tmp<-read.table(paste("PRS/scores_iPSYCHwoMINERvA_AllSNPs/iPSYCHwoMINERvA_AllSNPs.S", i, ".profile", sep = ""), header = TRUE)
	scores.all<-cbind(scores.all, scores.tmp[match(scores.tmp$IID, scores.all$IID),3:6])
	
}

colnames(scores.all)[seq(6,42,4)]<-paste("S", 1:10, sep = "")
colnames(scores.all)[seq(4,42,4)]<-paste("CNT_S", 1:10, sep = "")

### compare PRS across thresholds
pairs(scores.all[,seq(6,42,4)])
cor.all<-cor(scores.all[,seq(6,42,4)])
corrplot(cor.all, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, mar = c(5,5,5,5))


chart.Correlation(scores.all[,seq(6,42,4)], histogram=TRUE, pch=19)

pca<-read.table("01_Genotypes/Minerva.PCA.eigenvec", stringsAsFactors = FALSE)
pca<-pca[match(scores.all$IID, pca[,2]),]

### test pcas
scores.all$PHENO<-scores.all$PHENO-1

res.prs<-matrix(data = NA, nrow = 10, ncol = 2)
colnames(res.prs)<-c("Estimate", "P-value")
rownames(res.prs)<-paste("S", 1:10, sep = "")
pc1<-pca[,3]
pc2<-pca[,4]
pc3<-pca[,5]
pc4<-pca[,6]
pc5<-pca[,7]
for(i in 1:10){
	count<-scores.all[,seq(4,42,4)[i]]
	model<-glm(scores.all$PHENO ~ scores.all[,seq(6,42,4)[i]] + scores.all[,seq(4,42,4)[i]] + pca[,3] + pca[,4] + pca[,5] + pca[,6] + pca[,7],family=binomial(link='logit'))
	res.prs[i,1:2]<-summary(model)$coefficients["prs",c(1,4)]
}


barplot(-log10(res.prs[,2]), ylab = "-log10P", col = rainbow(10))

par(mfrow = c(2,5))
for(i in 1:10){
	d.case<-density(scores.all[which(scores.all$PHENO == 2),seq(6,42,4)[i]])
	d.con<-density(scores.all[which(scores.all$PHENO == 1),seq(6,42,4)[i]])
	plot(d.case, ylim = c(0,max(d.case$y, d.con$y)), xlab = "PRS", col = rgb(1,0,0,0.25), main = paste("S", i, sep = ""))
	polygon(d.case, col = rgb(1,0,0,0.25), border = rgb(1,0,0,0.25))
	polygon(d.con, col = rgb(0,1,0,0.25), border = rgb(0,1,0,0.25))
}
legend("topleft", c("ASD", "Control"), pch = 15, col = c(rgb(1,0,0,0.25), rgb(0,1,0,0.25)))

## Perform EWAS
load("QCed_Data.rda")

crosshyb<-read.csv("CrossHybridisingProbesPriceORWeksberg.csv", row.names = 1)
probes<-read.csv("/SNPsinProbesAnno.csv", row.names = 1)

## remove cross hybridising probes
remove<-match(crosshyb[,1], rownames(betas))
remove<-remove[which(is.na(remove) != TRUE)]
betas<-betas[-remove,]
## remove SNP probes 
probes<-probes[row.names(betas),]
betas<-betas[which(probes$Weksburg_CommonSNP_Af_within10bpSBE == "" | probes$Illumina_CommonSNP_Af_within10bpSBE == ""),]
betas<-betas[-grep("rs", rownames(betas)),]

scores.all<-scores.all[match(intersect(scores.all$IID, sampleSheet$Sample_Name), scores.all$IID),]
betas<-betas[, match(scores.all$IID, sampleSheet$Sample_Name)]
sampleSheet<-sampleSheet[match(scores.all$IID, sampleSheet$Sample_Name),]
pca<-pca[match(scores.all$IID, pca[,2]),]

## load covariates
counts<-read.csv("Analysis/HousemanEstimatedCellCounts.csv", stringsAsFactors = FALSE, row.names = 1)
counts<-counts[match(colnames(betas) ,rownames(counts)),]

scores<-read.csv("Analysis/QC/SmokingScores.csv", stringsAsFactors = FALSE, row.names = 1)
scores<-scores[match(colnames(betas), scores$Basename),]

res<-matrix(data = NA, nrow = nrow(betas), ncol = 30)
colnames(res)<-paste(paste("S", sort(rep(1:10,3)), sep = ""),c("P-value", "MeanDiff", "SE"), sep = ":")
rownames(res)<-rownames(betas)
for(i in 1:nrow(betas)){
	model<-lm(betas[i,] ~ scores.all$S1 + scores.all$CNT_S1 + pca[,3] + pca[,4] + pca[,5] + pca[,6] + pca[,7] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + pheno[, "fvagt"] + pheno[, "gest_week"])
	res[i,1:3]<-summary(model)$coefficients["scores.all$S1",c(1,2,4)]
	
	model<-lm(betas[i,] ~ scores.all$S2 + scores.all$CNT_S2 + pca[,3] + pca[,4] + pca[,5] + pca[,6] + pca[,7] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + pheno[, "fvagt"] + pheno[, "gest_week"])
	res[i,4:6]<-summary(model)$coefficients["scores.all$S2",c(1,2,4)]
	
	model<-lm(betas[i,] ~ scores.all$S3 + scores.all$CNT_S3 + pca[,3] + pca[,4] + pca[,5] + pca[,6] + pca[,7] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + pheno[, "fvagt"] + pheno[, "gest_week"])
	res[i,7:9]<-summary(model)$coefficients["scores.all$S3",c(1,2,4)]
	
	model<-lm(betas[i,] ~ scores.all$S4 + scores.all$CNT_S4 + pca[,3] + pca[,4] + pca[,5] + pca[,6] + pca[,7] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + pheno[, "fvagt"] + pheno[, "gest_week"])
	res[i,10:12]<-summary(model)$coefficients["scores.all$S4",c(1,2,4)]
	
	model<-lm(betas[i,] ~ scores.all$S5 + scores.all$CNT_S5 + pca[,3] + pca[,4] + pca[,5] + pca[,6] + pca[,7] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + pheno[, "fvagt"] + pheno[, "gest_week"])
	res[i,13:15]<-summary(model)$coefficients["scores.all$S5",c(1,2,4)]
	
	model<-lm(betas[i,] ~ scores.all$S6 + scores.all$CNT_S6 + pca[,3] + pca[,4] + pca[,5] + pca[,6] + pca[,7] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + pheno[, "fvagt"] + pheno[, "gest_week"])
	res[i,16:18]<-summary(model)$coefficients["scores.all$S6",c(1,2,4)]
	
	model<-lm(betas[i,] ~ scores.all$S7 + scores.all$CNT_S7 + pca[,3] + pca[,4] + pca[,5] + pca[,6] + pca[,7] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + pheno[, "fvagt"] + pheno[, "gest_week"])
	res[i,19:21]<-summary(model)$coefficients["scores.all$S7",c(1,2,4)]	
	
	model<-lm(betas[i,] ~ scores.all$S8 + scores.all$CNT_S8 + pca[,3] + pca[,4] + pca[,5] + pca[,6] + pca[,7] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + pheno[, "fvagt"] + pheno[, "gest_week"])
	res[i,22:24]<-summary(model)$coefficients["scores.all$S8",c(1,2,4)]
	
	model<-lm(betas[i,] ~ scores.all$S9 + scores.all$CNT_S9 + pca[,3] + pca[,4] + pca[,5] + pca[,6] + pca[,7] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + pheno[, "fvagt"] + pheno[, "gest_week"])
	res[i,25:27]<-summary(model)$coefficients["scores.all$S9",c(1,2,4)]
	
	model<-lm(betas[i,] ~ scores.all$S10 + scores.all$CNT_S10 + pca[,3] + pca[,4] + pca[,5] + pca[,6] + pca[,7] + as.factor(sampleSheet$gender.x.chr) + as.factor(sampleSheet$Sentrix_ID) + counts[,1] + counts[,2] + counts[,3] + counts[,4] + counts[,5] + counts[,6] + scores$scores_combined_A + pheno[, "fvagt"] + pheno[, "gest_week"])
	res[i,28:30]<-summary(model)$coefficients["scores.all$S10",c(1,2,4)]
}


library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
probeAnnot<-probeAnnot[rownames(betas),]
res<-cbind(res, probeAnnot)


write.csv(res, "Analysis/PRS_EWAS_CNT_5PCs_GA_BW_Sex_MaternalSmokingScoreCellComp.csv")


