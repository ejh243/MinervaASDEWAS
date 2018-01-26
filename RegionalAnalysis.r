##Perform sliding window regional analysis as described in Hannon et al. Genome Biology 2016.

calcCovMatrix<-function(x){
   ### x is a vector of correlations
   y<-vector(length = length(x))
   for(i in 1:length(x)){
      if(x[i] <= 1 & x[i] >= 0){
         y[i]<-x[i]*(3.25 + 0.75*x[i])
    } else {
      if(x[i] <= 0 & x[i] >= -0.5){
         y[i]<-x[i]*(3.27 + 0.71*x[i])
      }
    }
   
  }
  return(y)
}

brownsP<-function(covar, pval){
	## covar is vector of covariances between all pairs of tests NOTE not correlations
	## pval is vector of p values
	ntests<-length(pval)
	var_xsq<-sum(covar)*2 + 4*ntests 	# equation (3)
	exp_xsq<-2*ntests	#equation (2)

	## estimate parameters for chi square distribution using Brown's notation
	f = (2*(exp_xsq^2))/var_xsq
	c = var_xsq/(2*exp_xsq)

	##### NOTE: doesn't match Brown's answer but matches my answer by hand
	chi_sq<- -2*sum(log(pval))

	### to obtain p value
	test.stat<-chi_sq/c
	browns.pval<-1-pchisq(test.stat, df = f)
	return(c(browns.pval, test.stat, f))
}

## load DNA methylation data to calculate correlations
setwd("")
load("QCed_Data.rda")

crosshyb<-read.csv("450KAnno/CrossHybridisingProbesPriceORWeksberg.csv", row.names = 1)
probes<-read.csv("450KAnno/SNPsinProbesAnno.csv", row.names = 1)

## remove cross hybridising probes
remove<-match(crosshyb[,1], rownames(betas))
remove<-remove[which(is.na(remove) != TRUE)]
betas<-betas[-remove,]
## remove SNP probes as not a twin analysis
probes<-probes[row.names(betas),]
betas<-betas[which(probes$Weksburg_CommonSNP_Af_within10bpSBE == "" | probes$Illumina_CommonSNP_Af_within10bpSBE == ""),]

betas<-betas[-grep("rs", rownames(betas)),]

## load EWAS results
res<-read.csv("", row.names = 1, stringsAsFactors = FALSE)

res<-res[which(res$chr != "chrY"),]

betas<-betas[intersect(rownames(res), rownames(betas)),]
res<-res[rownames(betas),]
res<-res[order(res$pos),]


res$chr<-as.character(res$chr)
res$chr[which(res$chr == "chrX")]<-23
res$chr<-as.numeric(gsub("chr", "", res$chr))
res$UCSC_RefGene_Name<-as.character(res$UCSC_RefGene_Name)
windows<-c(1000,2000,5000)
out<-NULL
for(chr in 1:22){
	print(paste("CHR: ", chr))
	res.tmp<-res[which(res$chr == chr),]
	
	out.tmp<-matrix(data = NA, nrow = nrow(res.tmp), ncol = 2+length(windows)*7)
	colnames(out.tmp)<-c("Central Probe", "Chr", rep(c("Start Region", "End Region", "nProbes", "Gene Anno", "Min P", "Fishers P", "Browns P"), length(windows)))
	out.tmp[,1]<-rownames(res.tmp)	
	out.tmp[,2]<-chr
	for(i in 1:nrow(res.tmp)){

		chr<-res.tmp$chr[i]
		bp<-res.tmp$pos[i]

		for(j in 1:length(windows)){
		sub<-res.tmp[which(res.tmp$chr == chr & res.tmp$pos <= (bp+windows[j]) & res.tmp$pos >= (bp-windows[j])),]
		#sub<-sub[order(sub$pos),]
	
			out.tmp[i,(j*7)-4]<-min(sub$pos)
			out.tmp[i,(j*7)-3]<-max(sub$pos)
			out.tmp[i,(j*7)-2]<-nrow(sub)
			
			if(nrow(sub) > 1){
				out.tmp[i,(j*7)-1]<-paste(unique(unlist(strsplit(sub$UCSC_RefGene_Name, ";"))), collapse = ";")
				betas.sub<-betas[rownames(sub),]
				
				cor.mat<-c(unique(as.vector(cor(t(betas.sub)))))[-1]
				covar<-calcCovMatrix(cor.mat)

				pval<-sub$P.value
				out.tmp[i,(j*7)]<-min(pval)

				### compare to fishers p value
				fisher<-(sum(-log(pval))*2)
				out.tmp[i,(j*7)+1]<-1-pchisq(fisher, 2*length(pval))
				
				brownsp<-brownsP(covar, pval)
				out.tmp[i,(j*7)+2]<-brownsp[1]
			}
			
		}
	}
out<-rbind(out,out.tmp)
	
}


## save results for all tested regions
write.csv(out, "")

## filter results to unique regions tested
tmp2<-out[,c(1,2,10:16)]
tmp3<-out[,c(1,2,17:23)]
colnames(tmp2)<-colnames(out)[c(1:9)]
colnames(tmp3)<-colnames(out)[c(1:9)]
out2<-rbind(out[,c(1:9)], tmp2, tmp3)
out2<-unique(out2[,-1])
## filter regions to those with at least one DNA methylation site
out2<-out2[which(out2[,4] > 1),]
## sort by regional p value
out2<-out2[order(out2[,8]),]
