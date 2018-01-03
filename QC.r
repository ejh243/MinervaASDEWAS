## This script was used to load data from idat files into R, convert into beta values and perform a standard pipeline of data quality checks as outlined in the accompanying manuscript.
## Filepaths and sample IDs have been removed for security reasons, therfore it serves as a guide to the analysis

setwd("")
library(methylumi)
library(wateRmelon)

## load sample sheet and merge sentrix ids into unique identifier
sampleSheet<-read.csv("", stringsAsFactors = FALSE, skip = 8)
sampleSheet<-cbind(paste(sampleSheet$Sentrix_ID, sampleSheet$Sentrix_Position, sep = "_"), sampleSheet)
colnames(sampleSheet)[1]<-"Basename"
idatPath<-paste("idats/", sampleSheet$Sentrix_ID, "/", sampleSheet$Sentrix_ID, "_", sampleSheet$Sentrix_Position, sep ="")

mset450k <- methylumIDAT(sampleSheet$Basename[which(sampleSheet$Sentrix_ID == sampleSheet$Sentrix_ID[1])], idatPath=paste("idats/", sampleSheet$Sentrix_ID[1], "/", sep =""))
for(each in unique(sampleSheet$Sentrix_ID)[-c(1)]){
	mset450k.tmp <- methylumIDAT(sampleSheet$Basename[which(sampleSheet$Sentrix_ID == each)], idatPath=paste("idats/", each, "/", sep =""))
	mset450k<-combine(mset450k, mset450k.tmp)
}

## checking methylated and unmethylated signal intensities
m_intensities<-methylated(mset450k)
u_intensities<-unmethylated(mset450k)

M.median<-apply(m_intensities, 2, median)
U.median<-apply(u_intensities, 2, median)

sampleSheet<-sampleSheet[match(colnames(m_intensities), sampleSheet$Basename),]

## check profile of fully methylated and fully unmethylated control samples
pdf("Analysis/QC/Plots/MedianInstentiesPerSample_FullyMethylatedControls_Final.pdf")
par(mfrow = c(2,2))
boxplot(m_intensities[,which(sampleSheet$Sample_Group == "NegativCTRL")], names = sampleSheet$Sample_Plate[which(sampleSheet$Sample_Group == "NegativCTRL")], las = 2,  ylab = "M intensity", xlab = "", cex = 0.8)
abline(h = 3000, col = "red")
boxplot(u_intensities[,which(sampleSheet$Sample_Group == "NegativCTRL")],  names = sampleSheet$Sample_Plate[which(sampleSheet$Sample_Group == "NegativCTRL")], las = 2, ylab = "U intensity", xlab = "", cex = 0.8)
abline(h = 3000, col = "red")
boxplot(m_intensities[,which(sampleSheet$Sample_Group == "PositivCTRL")], names = sampleSheet$Sample_Plate[which(sampleSheet$Sample_Group == "PositivCTRL")], las = 2,  ylab = "M intensity", xlab = "", cex = 0.8)
abline(h = 3000, col = "red")
boxplot(u_intensities[,which(sampleSheet$Sample_Group == "PositivCTRL")],  names = sampleSheet$Sample_Plate[which(sampleSheet$Sample_Group == "PositivCTRL")], las = 2, ylab = "U intensity", xlab = "", cex = 0.8)
abline(h = 3000, col = "red")
dev.off()

## using the ten control probes to ensure the sodium bisulfite conversion was successful
greenChannel<-intensitiesByChannel(QCdata(mset450k))$Cy3
redChannel<-intensitiesByChannel(QCdata(mset450k))$Cy5

### these are fully methylated controls so closer to 1/100% the better the BS conversion.
### NOTE due to differnet chemistrys type I conversion controls and type II converstion controls need to be handled differently.
greenChannel<-greenChannel[grep("BS.Conversion", rownames(greenChannel)),]
redChannel<-redChannel[grep("BS.Conversion", rownames(redChannel)),]
BScon1<-rbind(greenChannel[c("BS.Conversion.I.C1", "BS.Conversion.I.C2", "BS.Conversion.I.C3"),], redChannel[c("BS.Conversion.I.C4", "BS.Conversion.I.C5", "BS.Conversion.I.C6"),])
BScon2<-redChannel[c("BS.Conversion.II.1", "BS.Conversion.II.2", "BS.Conversion.II.3", "BS.Conversion.II.4"),]

### for type II probes red channel is M and green channel is U so easy conversion into beta values
BScon2.betas<-redChannel[c("BS.Conversion.II.1", "BS.Conversion.II.2", "BS.Conversion.II.3", "BS.Conversion.II.4"),]/(redChannel[c("BS.Conversion.II.1", "BS.Conversion.II.2", "BS.Conversion.II.3", "BS.Conversion.II.4"),] + greenChannel[c("BS.Conversion.II.1", "BS.Conversion.II.2", "BS.Conversion.II.3", "BS.Conversion.II.4"),])

### for type I both methylated and unmethylated in the same channel
BScon1.betas<-rbind(greenChannel[c("BS.Conversion.I.C1", "BS.Conversion.I.C2", "BS.Conversion.I.C3"),], redChannel[c("BS.Conversion.I.C4", "BS.Conversion.I.C5", "BS.Conversion.I.C6"),])/(rbind(greenChannel[c("BS.Conversion.I.C1", "BS.Conversion.I.C2", "BS.Conversion.I.C3"),], redChannel[c("BS.Conversion.I.C4", "BS.Conversion.I.C5", "BS.Conversion.I.C6"),]) + rbind(greenChannel[c("BS.Conversion.I.U1", "BS.Conversion.I.U2", "BS.Conversion.I.U3"),], redChannel[c("BS.Conversion.I.U4", "BS.Conversion.I.U5", "BS.Conversion.I.U6"),]))

BScon.betas<-rbind(BScon1.betas, BScon2.betas)
BScon.beta.median<-apply(BScon.betas, 2, median)*100

pdf("Analysis/QC/Plots/Histogram_MedianBSConversion_Final.pdf", width = 18, height = 8)
hist(BScon.beta.median, xlab = "Median % BS conversion", main = "")
abline(v = 80, col = "red")
dev.off()

pdf("Analysis/QC/Plots/Boxplot_MedianBSConversion_Final.pdf", width = 18, height = 8)
boxplot(BScon.beta.median ~ sampleSheet$Sample_Plate, ylab = "Median % BS conversion", col = rainbow(14), las = 2)
abline(h = 80, col = "red")
dev.off()

### compare sample intensity with bisulfite conversion
col_points<-rainbow(14)[factor(sampleSheet$Sample_Plate)]
pch_points<-rep(16, length(col_points))
pch_points[which(sampleSheet$Sample_Group == "Negativ" | sampleSheet$Sample_Group == "Positiv")]<-17
pdf("Analysis/QC/Plots/Scatterplot_MedianBSConversion_Intensities_Final.pdf", width = 18, height = 8)
par(mfrow = c(1,2))
plot(M.median, BScon.beta.median, xlab = "Median M intensity", ylab = "Median % BS conversion", col = col_points, pch = pch_points)
plot(U.median, BScon.beta.median,  xlab = "Median U intensity", ylab = "Median % BS conversion", col = col_points, pch = pch_points)
legend("bottomright", levels(factor(sampleSheet$Sample_Plate)), col = rainbow(14), pch = 16)
dev.off()

sampleSheet<-cbind(sampleSheet, BScon.beta.median)

### remove control samples & those failing bisulfite conversion
mset450k<-mset450k[,-which(sampleSheet$Sample_Group == "NegativCTRL")]
sampleSheet<-sampleSheet[-which(sampleSheet$Sample_Group == "NegativCTRL"),]
mset450k<-mset450k[,-which(sampleSheet$Sample_Group == "PositivCTRL")]
sampleSheet<-sampleSheet[-which(sampleSheet$Sample_Group == "PositivCTRL"),]

### remove those failing bisulfite conversion
mset450k<-mset450k[,-which(sampleSheet$BScon.beta.median < 80)]
sampleSheet<-sampleSheet[-which(sampleSheet$BScon.beta.median < 80),]


### check for duplicate samples
meth<-betas(mset450k)
meth.rs<-meth[grep("rs", rownames(meth)),]
nNA<-colSums(is.na(meth.rs))
pdf("Analysis/QC/Plots/Boxplot_nNAinSNPProbes_Final.pdf")
boxplot(nNA ~ sampleSheet$Sample_Plate, col = rainbow(14), ylab = "N missing values", las = 2)
dev.off()

cor.mat<-cor(meth.rs, use = "pairwise.complete.obs")
for(i in 1:ncol(cor.mat)){

	cor.mat[i,i]<-NA
}

### ignore anything with more than 5 missing data points

cor.mat[which(nNA > 5),]<-NA
cor.mat[,which(nNA > 5)]<-NA

maxCor<-apply(cor.mat, 1, max, na.rm = TRUE)

pdf("Analysis/QC/Plots/SNPProbeCorrelations_Final.pdf")
hist(maxCor, xlab = "Maximum SNP correlation with all other samples")
dev.off()

##### age calculator
library(WGCNA)
library(sqldf)
library(impute)
library(RPMM)
#install the Bioconductor installer
#install.packages("BiocInstaller",repos="http://www.bioconductor.org/packages/2.13/bioc")
#install “impute” from Bioconductor
#source("http://bioconductor.org/biocLite.R")
#biocLite("impute")
## need unnormalised data as script includes normalisation
source("/data/Minerva/FromEilis/NORMALIZATION.R")

#Age transformation and probe annotation functions
trafo= function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo= function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
probeAnnotation21kdatMethUsed=read.csv("/data/Minerva/FromEilis/probeAnnotation21kdatMethUsed.csv")
probeAnnotation27k=read.csv("/data/Minerva/FromEilis/datMiniAnnotation27k.csv")
datClock=read.csv("/data/Minerva/FromEilis/AdditionalFile3.csv")

#Read in the DNA methylation data (beta values)
# For a small file, e.g. measured on the 27k platform you could just use read.csv. 
# But for large files, e.g. those measured on the 450K platform, I recommend you use read.csv.sql.

#dat0=read.csv.sql("/data/Minerva/FromEilis/MethylationDataExample55.csv") ; ### load in data

full.data.betas<-betas(mset450k)
full.data.betas<-as.data.frame(full.data.betas)
dat0<-as.data.frame("ProbeIDs" = rownames(full.data.betas), full.data.betas) ## first column is list of probe names

full.data.betas<-betas(mset450k)
full.data.betas<-as.data.frame(full.data.betas)
dat0<-cbind(rownames(full.data.betas), full.data.betas) ## first column is list of probe names
colnames(dat0)[1]<-"ProbeIDs"
dat0[,1]<-as.character(dat0[,1])
nSamples=dim(dat0)[[2]]-1
nProbes= dim(dat0)[[1]]
# the following command may not be needed. But it is sometimes useful when you use read.csv.sql
#dat0[,1]= gsub(x=dat0 [,1],pattern="\"",replacement="") 
#Create a log file which will be output into your directory
# The code looks a bit complicated because it serves to create a log file (for error checks etc).
# It will automatically create a log file.
file.remove("Analysis/QC/LogFile.txt")
file.create("Analysis/QC/LogFile.txt")
DoNotProceed=FALSE
cat(paste( "The methylation data set contains", nSamples, "samples (e.g. arrays) and ", nProbes, " probes."),file="Analysis/QC/LogFile.txt")
if (nSamples==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be no samples.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql . Samples correspond to columns in that file  ."), file="Analysis/QC/LogFile.txt",append=TRUE) } 
if (nProbes==0) {DoNotProceed=TRUE; cat(paste( "\n ERROR: There must be a data input error since there seem to be zero probes.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql  CpGs correspond to rows.")   , file="Analysis/QC/LogFile.txt",append=TRUE) } 
if (  nSamples > nProbes  ) { cat(paste( "\n MAJOR WARNING: It worries me a lot that there are more samples than CpG probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose the data and then resubmit them? In any event, I will proceed with the analysis."),file="Analysis/QC/LogFile.txt",append=TRUE) }
if (  is.numeric(dat0[,1]) ) { DoNotProceed=TRUE; cat(paste( "\n Error: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="Analysis/QC/LogFile.txt",append=TRUE)  } 
if (  !is.character(dat0[,1]) ) {  cat(paste( "\n Major Warning: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]  ),file="Analysis/QC/LogFile.txt",append=TRUE)  } 
datout=data.frame(Error=c("Input error. Please check the log file for details","Please read the instructions carefully."), Comment=c("", "email Steve Horvath."))
if ( ! DoNotProceed ) {
nonNumericColumn=rep(FALSE, dim(dat0)[[2]]-1)
for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
if (  sum(nonNumericColumn) >0 ) { cat(paste( "\n MAJOR WARNING: Possible input error. The following samples contain non-numeric beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure this makes sense.\n" ),file="Analysis/QC/LogFile.txt",append=TRUE)  } 
XchromosomalCpGs=as.character(probeAnnotation27k$Name[probeAnnotation27k$Chr=="X"])
selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
selectXchromosome[is.na(selectXchromosome)]=FALSE
meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
if (   sum(selectXchromosome) >=500 )  {
meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE)) }
if (  sum(is.na(meanXchromosome)) >0 ) { cat(paste( "\n \n Comment: There are lots of missing values for X chromosomal probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these samples.\n " ),file="Analysis/QC/LogFile.txt",append=TRUE)  } 

match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
if  ( sum( is.na(match1))>0 ) { 
missingProbes= probeAnnotation21kdatMethUsed$Name[!is.element( probeAnnotation21kdatMethUsed$Name , dat0[,1])]    
DoNotProceed=TRUE; cat(paste( "\n \n Input error: You forgot to include the following ", length(missingProbes), " CpG probes (or probe names):\n ", paste( missingProbes, sep="",collapse=", ")),file="Analysis/QC/LogFile.txt",append=TRUE)  } 

#STEP 2: Restrict the data to 21k probes and ensure they are numeric
match1=match(probeAnnotation21kdatMethUsed$Name , dat0[,1])
if  ( sum( is.na(match1))>0 ) stop(paste(sum( is.na(match1)), "CpG probes cannot be matched"))
dat1= dat0[match1,]
asnumeric1=function(x) {as.numeric(as.character(x))}
dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)

#STEP 3: Create the output file called datout
set.seed(1)
# Do you want to normalize the data (recommended)?
normalizeData=TRUE
source("/data/Minerva/FromEilis/StepwiseAnalysis.txt")
# STEP 4: Output the results 
if (  sum(  datout$Comment  != "" )   ==0 ) { cat(paste( "\n The individual samples appear to be fine. "),file="Analysis/QC/LogFile.txt",append=TRUE)  } 
if (  sum(  datout$Comment != "" )   >0 ) { cat(paste( "\n Warnings were generated for the following samples.\n", datout[,1][datout$Comment != ""], "\n Hint: Check the output file for more details."),file="Analysis/QC/LogFile.txt",append=TRUE)  } 
} 
# output the results into the directory
write.table(datout,"Analysis/QC/AgeCalculator_Output.csv", row.names=F, sep="," )

### merge with phenotype info
datout<-cbind(datout, sampleSheet[match(rownames(datout), sampleSheet$Basename),])
write.table(datout,"Analysis/QC/AgeCalculator_Output_Final.csv", row.names=F, sep="," )

## load and match additional phenotype information
pheno<-read.csv("", stringsAsFactors = FALSE)
pheno.cases<-pheno[,c(21, 18, grep("case", colnames(pheno)))]
pheno.cont<-pheno[,c(22, 18, grep("cont", colnames(pheno)))]
colnames(pheno.cases)<-c("participant_id", gsub("_case", "", colnames(pheno.cases))[-1])
colnames(pheno.cont)<-c("participant_id", gsub("_cont", "", colnames(pheno.cont))[-1])
pheno.cont<-pheno.cont[,colnames(pheno.cases)]

pheno<-rbind(pheno.cases, pheno.cont)
pheno<-pheno[match(sampleSheet$Sample_Name, gsub("MMXII_iPSYCH_", "", pheno$participant_id)),]

### compare reported sex with predicted sex
table(datout$Gender, datout$predictedGender)

pdf("Analysis/QC/Plots/Histogram_DNAmAge.pdf")
hist(datout$DNAmAge, xlab = "Predicted methylation age")
boxplot(datout$DNAmAge ~ datout$Sample_Plate, ylab = "Predicted methylation age", col = rainbow(14), las = 2)
dev.off()


### check gender using multidimensional scaling of data from probes on the X and Y chromosomes
betas<-betas(mset450k)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
probeAnnot<-probeAnnot[rownames(betas),]

### X chr
probeAnnot.x<-probeAnnot[which(probeAnnot[,1] == "chrX"),]
betas.x<-betas[rownames(probeAnnot.x),]

d<-dist(t(betas.x))
fit<-cmdscale(d, k = 2)

pdf("Analysis/QC/Plots/MDS_SexChromosomes_GenderChecks_Final.pdf", paper = "a4r", width = 0, height = 0)
par(mfrow = c(1,2))
sex_palette<-cbind(c("M", "F"), c("blue", "magenta"))
sex_col<-sex_palette[match(sampleSheet$Gender, sex_palette[,1]),2]
plot(fit[,1], fit[,2], xlab = "MDS Co-ordinate 1", ylab = "MDS Co-ordinate 2", pch = 18, col = sex_col, main = "X Chromosome")
legend("topright", horiz = TRUE, c("M", "F"), col = c("blue", "magenta"), pch = 18)

gender.x.chr<-rep("Unsure", length = nrow(sampleSheet))
gender.x.chr[which(fit[,1] > 0)]<-"F"
gender.x.chr[which(fit[,1] < 0)]<-"M"

fit.x<-fit
### Y chr
probeAnnot.y<-probeAnnot[which(probeAnnot[,1] == "chrY"),]
betas.y<-betas[rownames(probeAnnot.y),]

d<-dist(t(betas.y))
fit<-cmdscale(d, k = 2)

sex_palette<-cbind(c("M", "F"), c("blue", "magenta"))
sex_col<-sex_palette[match(sampleSheet$Gender, sex_palette[,1]),2]
plot(fit[,1], fit[,2], xlab = "MDS Co-ordinate 1", ylab = "MDS Co-ordinate 2", pch = 18, col = sex_col, main = "Y Chromosome")
legend("topright", horiz = TRUE, c("M", "F"), col = c("blue", "magenta"), pch = 18)
dev.off()

gender.y.chr<-vector(length = nrow(sampleSheet))
gender.y.chr[which(fit[,1] > 0)]<-"F"
gender.y.chr[which(fit[,1] < 0)]<-"M"

sampleSheet<-cbind(sampleSheet,gender.x.chr, gender.y.chr)

pdf("Analysis/QC/Plots/Scatterplot_MDS_xchr_ychr.pdf", height = 8, width = 8)
plot(fit.x[,1], fit[,1], pch = 16, xlab = "X chr", ylab = "Y chr", col = c("blue", "magenta")[factor(datout$Gender)])
abline(v = 0)
abline(h = 0)
dev.off()

### identify which samples do not match gender prediction
which(gender.x.chr != gender.y.chr)
mset450k<-mset450k[,which(gender.x.chr == gender.y.chr)]
sampleSheet<-sampleSheet[which(gender.x.chr == gender.y.chr),]


###pfilter
mset.pf<-pfilter(mset450k)

###normalisation
mset.dasen<-dasen(mset.pf)

betas<-betas(mset.dasen)

### exclude two samples  subsquently identified as swaps with genotype data
exclude<-c("", "")
betas<-betas[,-match(exclude, sampleSheet$Sample_Name)]
sampleSheet<-sampleSheet[-match(exclude, sampleSheet$Sample_Name),]
save(betas, sampleSheet, datout, file = "QCed_Data.rda")
