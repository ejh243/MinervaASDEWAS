## This script was used to perform the meta-analysis combining the p values from EWAS performed in three different cohorts

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(qqman)

setwd("")
load("") ## load seed.SS object
load("") ## load simons.SS object
minerva<-read.csv("Analysis/CaseControl_Sex_Batch_HousemanCellComp_SmokingScore_BY_BM_GA_Urbanicity.csv", stringsAsFactors = FALSE, row.names = 1)

probes<-table(c(rownames(minerva), rownames(seed.SS), rownames(simons.SS)))
table(probes)

minerva<-minerva[names(probes)[which(probes > 1)],]
seed.SS<-seed.SS[names(probes)[which(probes > 1)],]
simons.SS<-simons.SS[names(probes)[which(probes > 1)],]

### look at overlap of top probes

par(mfrow = c(3,2))
plot(minerva$MeanDiff[which(minerva$P.value < 5e-5)], seed.SS$MeanDiff[which(minerva$P.value < 5e-5)], xlab = "Minerva", ylab = "seed", main = paste(length(which(sign(minerva$MeanDiff[which(minerva$P.value < 5e-5)]) == sign(seed.SS$MeanDiff[which(minerva$P.value < 5e-5)]))), "/", length(which(minerva$P.value < 5e-5)), " P = ", signif(binom.test(length(which(sign(minerva$MeanDiff[which(minerva$P.value < 5e-5)]) == sign(seed.SS$MeanDiff[which(minerva$P.value < 5e-5)]))), length(which(minerva$P.value < 5e-5)))$p.value,3), sep = ""), pch = 16)
abline(v = 0)
abline(h = 0)
plot(minerva$MeanDiff[which(minerva$P.value < 5e-5)], simons.SS$MeanDiff[which(minerva$P.value < 5e-5)], xlab = "Minerva", ylab = "simons", main = paste(length(which(sign(minerva$MeanDiff[which(minerva$P.value < 5e-5)]) == sign(simons.SS$MeanDiff[which(minerva$P.value < 5e-5)]))), "/", length(which(minerva$P.value < 5e-5)), " P = ", signif(binom.test(length(which(sign(minerva$MeanDiff[which(minerva$P.value < 5e-5)]) == sign(simons.SS$MeanDiff[which(minerva$P.value < 5e-5)]))), length(which(minerva$P.value < 5e-5)))$p.value,3), sep = ""), pch = 16)
abline(v = 0)
abline(h = 0)

plot(seed.SS$MeanDiff[which(seed.SS$pvalue < 5e-5)], minerva$MeanDiff[which(seed.SS$pvalue < 5e-5)], xlab = "seed", ylab = "Minerva", main = paste(length(which(sign(seed.SS$MeanDiff[which(seed.SS$pvalue < 5e-5)]) == sign(minerva$MeanDiff[which(seed.SS$pvalue < 5e-5)]))), "/", length(which(seed.SS$pvalue < 5e-5)), " P = ", signif(binom.test(length(which(sign(seed.SS$MeanDiff[which(seed.SS$pvalue < 5e-5)]) == sign(minerva$MeanDiff[which(seed.SS$pvalue < 5e-5)]))), length(which(seed.SS$pvalue < 5e-5)))$p.value,3), sep = ""), pch = 16)
abline(v = 0)
abline(h = 0)
plot(seed.SS$MeanDiff[which(seed.SS$pvalue < 5e-5)], simons.SS$MeanDiff[which(seed.SS$pvalue < 5e-5)], xlab = "seed", ylab = "simons", main = paste(length(which(sign(seed.SS$MeanDiff[which(seed.SS$pvalue < 5e-5)]) == sign(simons.SS$MeanDiff[which(seed.SS$pvalue < 5e-5)]))), "/", length(which(seed.SS$pvalue < 5e-5)), " P = ", signif(binom.test(length(which(sign(seed.SS$MeanDiff[which(seed.SS$pvalue < 5e-5)]) == sign(simons.SS$MeanDiff[which(seed.SS$pvalue < 5e-5)]))), length(which(seed.SS$pvalue < 5e-5)))$p.value,3), sep = ""), pch = 16)
abline(v = 0)
abline(h = 0)

plot(simons.SS$MeanDiff[which(simons.SS$pvalue < 5e-5)], minerva$MeanDiff[which(simons.SS$pvalue < 5e-5)], xlab = "simons", ylab = "Minerva", main = paste(length(which(sign(simons.SS$MeanDiff[which(simons.SS$pvalue < 5e-5)]) == sign(minerva$MeanDiff[which(simons.SS$pvalue < 5e-5)]))), "/", length(which(simons.SS$pvalue < 5e-5)), " P = ", signif(binom.test(length(which(sign(simons.SS$MeanDiff[which(simons.SS$pvalue < 5e-5)]) == sign(minerva$MeanDiff[which(simons.SS$pvalue < 5e-5)]))), length(which(simons.SS$pvalue < 5e-5)))$p.value,3), sep = ""), pch = 16)
abline(v = 0)
abline(h = 0)
plot(simons.SS$MeanDiff[which(simons.SS$pvalue < 5e-5)], seed.SS$MeanDiff[which(simons.SS$pvalue < 5e-5)], xlab = "simons", ylab = "seed", main = paste(length(which(sign(simons.SS$MeanDiff[which(simons.SS$pvalue < 5e-5)]) == sign(seed.SS$MeanDiff[which(simons.SS$pvalue < 5e-5)]))), "/", length(which(simons.SS$pvalue < 5e-5)), " P = ", signif(binom.test(length(which(sign(simons.SS$MeanDiff[which(simons.SS$pvalue < 5e-5)]) == sign(seed.SS$MeanDiff[which(simons.SS$pvalue < 5e-5)]))), length(which(simons.SS$pvalue < 5e-5)))$p.value,3), sep = ""), pch = 16)
abline(v = 0)
abline(h = 0)



par(mfrow = c(3,2))
plot(-log10(minerva$P.value[which(minerva$P.value < 5e-5)]), -log10(seed.SS$pvalue[which(minerva$P.value < 5e-5)]), xlab = "Minerva", ylab = "seed", pch = 16)
abline(v = 0)
abline(h = 0)
plot(-log10(minerva$P.value[which(minerva$P.value < 5e-5)]), -log10(simons.SS$pvalue[which(minerva$P.value < 5e-5)]), xlab = "Minerva", ylab = "simons", pch = 16)
abline(v = 0)
abline(h = 0)


### use fisher's method to combine p-values

res<-cbind(minerva[match(names(probes)[which(probes == 3)], rownames(minerva)), c("MeanDiff", "SE", "P.value")], seed.SS[match(names(probes)[which(probes == 3)], rownames(seed.SS)) , c("MeanDiff", "SE", "pvalue")], simons.SS[match(names(probes)[which(probes == 3)], rownames(simons.SS)), c("MeanDiff", "SE", "pvalue")])
res<-cbind(res, NA, 3)
colnames(res)<-c("Minerva:MeanDiff", "Minerva:SE", "Minerva:P", "seed:MeanDiff", "seed:SE", "seed:P","simons:MeanDiff", "simons:SE", "simons:P","Fishers_P", "n_studies")
rownames(res)<-names(probes)[which(probes ==3)]

### start with probes annotated in all 3 datasets
for(each in names(probes)[which(probes == 3)]){

	pval<-res[each,c(3,6,9)]

	### compare to fishers p value
	fisher<-(sum(-log(pval))*2)
	fisherp<-1-pchisq(fisher, 2*length(pval))
	
	res[each,10]<-fisherp
}

write.csv(res, "Analysis/CaseControl_ASD_MetaAnalysis_FishersP.csv")

### then consider probes tested in 2 datasets

res.2<-cbind(minerva[match(names(probes)[which(probes == 2)], rownames(minerva)), c("MeanDiff", "SE", "P.value")], seed.SS[match(names(probes)[which(probes == 2)], rownames(seed.SS)) , c("MeanDiff", "SE", "pvalue")], simons.SS[match(names(probes)[which(probes == 2)], rownames(simons.SS)), c("MeanDiff", "SE", "pvalue")])
res.2<-cbind(res.2, NA, 2)
colnames(res.2)<-c("Minerva:MeanDiff", "Minerva:SE", "Minerva:P", "seed:MeanDiff", "seed:SE", "seed:P","simons:MeanDiff", "simons:SE", "simons:P","Fishers_P", "n_studies")
rownames(res.2)<-names(probes)[which(probes == 2)]
res<-rbind(res,res.2)

for(each in names(probes)[which(probes == 2)]){
	if(is.na(minerva$MeanDiff[match(each, rownames(minerva))])){
	pval<-res[each,c(6,9)]
	} else {
		if(is.na(seed.SS$MeanDiff[match(each, rownames(seed.SS))])){		
			pval<-res[each,c(3,9)]
		} else {
			if(is.na(simons.SS$MeanDiff[match(each, rownames(simons.SS))])){
				
				pval<-res[each,c(3,6)]
			}
		}	
	}
	fisher<-(sum(-log(pval))*2)
	fisherp<-1-pchisq(fisher, 2*length(pval))
	res[each,10]<-fisherp
}


write.csv(res, "Analysis/CaseControl_ASD_MetaAnalysis_FishersP.csv")


data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Locations)
probeAnnot<-probeAnnot[rownames(res),]
res<-cbind(res, probeAnnot)
probeAnnot<-as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other)
probeAnnot<-probeAnnot[rownames(res),]
res<-cbind(res, probeAnnot)
write.csv(res, "Analysis/CaseControl_ASD_MetaAnalysis_FishersP.csv")

qq(res$Fishers_P)
res$chr<-gsub("chr", "", res$chr)
res$chr[which(res$chr == "X")]<-23
res$chr<-as.numeric(res$chr)

##For manhattan plot redefine function to use symbools to indicate when all same direction of effects
manhattan<-function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10",
    "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05),
    genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE,
    annotatePval = NULL, annotateTop = TRUE, pch = 20, ...)
{
    CHR = BP = P = index = NULL
    if (!(chr %in% names(x)))
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x)))
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x)))
        stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x)))
        warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!is.numeric(x[[chr]]))
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]]))
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]]))
        stop(paste(p, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
    if (!is.null(x[[snp]]))
        d = transform(d, SNP = x[[snp]])
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
    if (logp) {
        d$logp <- -log10(d$P)
    }
    else {
        d$logp <- d$P
    }
    d$pos = NA
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        d$pos = d$BP
        ticks = floor(length(d$pos))/2 + 1
        xlabel = paste("Chromosome", unique(d$CHR), "position")
        labs = ticks
    }
    else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            }
            else {
                lastbase = lastbase + tail(subset(d, index ==
                  i - 1)$BP, 1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP +
                  lastbase
            }
            ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index ==
                i, ]$pos))/2 + 1)
        }
        xlabel = "Chromosome"
        labs <- unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i",
        las = 1, pch = pch, xlim = c(xmin, xmax), ylim = c(0,
            ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
    dotargs <- list(...)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in%
        names(dotargs)]))
    if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
            if (length(chrlabs) == length(labs)) {
                labs <- chrlabs
            }
            else {
                warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
            }
        }
        else {
            warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
    }
    if (nchr == 1) {
        axis(1, ...)
    }
    else {
        axis(1, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logp, pch = pch, col = col[1], ...))
    }
    else {
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], points(pos,
                logp, col = col[icol], pch = pch, ...))
            icol = icol + 1
        }
    }
    if (suggestiveline)
        abline(h = suggestiveline, col = "blue")
    if (genomewideline)
        abline(h = genomewideline, col = "red")
    if (!is.null(highlight)) {
        if (any(!(highlight %in% d$SNP)))
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight = d[which(d$SNP %in% highlight), ]
        with(d.highlight, points(pos, logp, col = "green3", pch = pch,
            ...))
    }
    if (!is.null(annotatePval)) {
        topHits = subset(d, P <= annotatePval)
        par(xpd = TRUE)
        if (annotateTop == FALSE) {
            with(subset(d, P <= annotatePval), textxy(pos, -log10(P),
                offset = 0.625, labs = topHits$SNP, cex = 0.45),
                ...)
        }
        else {
            topHits <- topHits[order(topHits$P), ]
            topSNPs <- NULL
            for (i in unique(topHits$CHR)) {
                chrSNPs <- topHits[topHits$CHR == i, ]
                topSNPs <- rbind(topSNPs, chrSNPs[1, ])
            }
            textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625,
                labs = topSNPs$SNP, cex = 0.5, ...)
        }
    }
    par(xpd = FALSE)
}

pchInd<-rep(2, nrow(res))
pchInd[which(abs(rowSums(sign(res[,c(1,4,7)]), na.rm = TRUE)) == rowSums(!is.na(res[,c(1,4,7)])))]<-16

manhattan(res, p = "Fishers_P", chr = "chr", bp = "pos", suggestiveline = -log10(5e-5), genomewideline = -log10(1e-7), pch = pchInd)
