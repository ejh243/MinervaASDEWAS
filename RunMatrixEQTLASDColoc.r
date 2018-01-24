library(MatrixEQTL)


setwd("")
base.dir = find.package("MatrixEQTL")
useModel = modelLINEAR;
pvOutputThreshold = 1;
cisDist<-500000
errorCovariance = numeric();

## only first three ar genome-wide significant
regions<-read.table("Analysis/AdditionalAutismRegions.txt", header = TRUE, row.names = 1)
regions<-regions[1:3,]

covariates_file_name = "MatrixEQTL/Covariates_Sex_5PCs.txt"; 

##load covariate data
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile(covariates_file_name);


for(chr in unique(regions$CHR)){
	SNP_file_name = paste("MatrixEQTL/Genotypes_MinGenoCount5_Chr", chr, ".txt", sep = "")
	snps_location_file_name = paste("MatrixEQTL/Genotypes_MapInfo_Chr", chr, ".txt", sep = "")
	
	 ## load SNP data
	snps = SlicedData$new();
	snps$fileDelimiter = "\t";      # the TAB character
	snps$fileOmitCharacters = "NA"; # denote missing values;
	snps$fileSkipRows = 1;          # one row of column labels
	snps$fileSkipColumns = 1;       # one column of row labels
	snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
	snps$LoadFile( SNP_file_name );

	snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
	snpspos[,1]<-rownames(snps)

for(i in which(regions$CHR == chr)){

	start<-regions$BP[i]-250000
	stop<-regions$BP[i]+250000
	expression_file_name<-paste("MatrixEQTL/Methylation_ProbesForASDColocAnalysis_Chr", chr, "_", start, "_", stop, ".txt", sep = "")

	#gene_location_file_name = paste(base.dir, "/data/geneloc.txt", sep="");
	## if no covariates set to  character()
	output_file_name = paste("MatrixEQTL/Output/AllmQTLsforASDColocAnalysis_mQTL_chr", chr,"_", start, "_", stop, ".txt", sep = "")

	con <- file(expression_file_name) 
		
		if (length(readLines(con)) > 1){
			## load methylation data
			gene = SlicedData$new();
			gene$fileDelimiter = "\t";      # the TAB character
			gene$fileOmitCharacters = "NA"; # denote missing values;
			gene$fileSkipRows = 1;          # one row of column labels
			gene$fileSkipColumns = 1;       # one column of row labels
			gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
			gene$LoadFile( expression_file_name);

		library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
		data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
		probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
		probeAnnot<-cbind(rownames(gene), probeAnnot[rownames(gene), c("chr", "pos", "pos")])
		probeAnnot[,4]<-probeAnnot[,4]+1
		genepos<-probeAnnot

			### run eqtls
			me = Matrix_eQTL_main(
				snps = snps,
				gene = gene,
				cvrt = cvrt,
				output_file_name = output_file_name,
				pvOutputThreshold = 1,
				output_file_name.cis = output_file_name,
				pvOutputThreshold.cis = pvOutputThreshold,
				snpspos = snpspos, 
				genepos = genepos,
				cisDist = cisDist,
				useModel = useModel, 
				errorCovariance = errorCovariance, 
				verbose = TRUE,
				pvalue.hist = "qqplot",
				min.pv.by.genesnp = FALSE,
				noFDRsaveMemory = FALSE)
				
		}
	}
}
