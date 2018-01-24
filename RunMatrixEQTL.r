## Run matrixQTL to perform whole genome scan for DNA methylation quantitative trait loci. 
## Analysis is performed by analysising all genetic variants on each chromosome in turn against all DNA methylation sites. 

library(MatrixEQTL)

setwd("")
base.dir = find.package("MatrixEQTL")
useModel = modelLINEAR; 

covariates_file_name = "MatrixEQTL/Covariates_Sex_5PCs.txt"; 
expression_file_name = "MatrixEQTL/Methylation.txt"

## threshold of results to save
pvOutputThreshold = 1e-8;
errorCovariance = numeric();

##load covariate data
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
cvrt$LoadFile(covariates_file_name);

## load methylation data
gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
gene$LoadFile( expression_file_name);

for(chr in 1:22){
	SNP_file_name = paste("MatrixEQTL/Genotypes_MinGenoCount5_Chr", chr, ".txt", sep = "")
	output_file_name = paste("MatrixEQTL/Output/Minerva_mQTL_chr", chr, ".txt", sep = "")

	 ## load SNP data
	snps = SlicedData$new();
	snps$fileDelimiter = "\t";      # the TAB character
	snps$fileOmitCharacters = "NA"; # denote missing values;
	snps$fileSkipRows = 1;          # one row of column labels
	snps$fileSkipColumns = 1;       # one column of row labels
	snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
	snps$LoadFile( SNP_file_name );

	### run eqtls
	me = Matrix_eQTL_main(
		snps = snps,
		gene = gene,
		cvrt = cvrt,
		output_file_name = output_file_name,
		pvOutputThreshold = pvOutputThreshold,
		useModel = useModel, 
		errorCovariance = errorCovariance, 
		verbose = TRUE,
		pvalue.hist = "qqplot",
		min.pv.by.genesnp = FALSE,
		noFDRsaveMemory = FALSE)
		
}
