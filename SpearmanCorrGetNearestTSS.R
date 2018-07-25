### R code to correlate gene expression and DNA methylation
## By using getNearestTSS from FDb.InfiniumMethylation.hg19
load("../GeneExp.Rda")## RSEM normalized read counts
load("../BMIQ_Meth.Rda")## BMIQ normalized Beta value
load("../../probe.features.Rda")## Illumina gene name and short annotation
#SampleId <- intersect(colnames(geneExp), colnames(BMIQ.Meth))
geneExp <- round(geneExp);geneExp <- log2(geneExp+1)
rownames(geneExp) <- gsub(pattern = "SLC35E2\\|728661", "SLC35E2B\\|728661", x = rownames(geneExp))
geneExp <- geneExp[-grep("\\?", rownames(geneExp)),] ## Remove genes which HGNC name is known i.e. ?
geneExp$gene <- sapply(strsplit(rownames(geneExp),"\\|"),'[[',1)
rownames(geneExp) <- geneExp$gene
#mergeExp <- merge(geneExp, DiffExp, by="gene")
#rownames(mergeExp) <- mergeExp$gene
mergeMeth <- merge(BMIQ.Meth, probe.features, by="row.names")
rownames(mergeMeth) <- mergeMeth$Row.names
#####################################################################
library(FDb.InfiniumMethylation.hg19)
hm450 <- get450k()
df1 <- data.frame(seqnames=names(hm450),
                  starts=start(hm450),
                  ends=end(hm450),
                  names=c(rep(".", length(hm450))),
                  strand=strand(hm450))
probenames <- df1$seqnames
probes <- hm450[probenames]
TSS.probes <- getNearestTSS(probes)
TSS.probes$queryHits <- NULL; TSS.probes$subjectHits <- NULL
TSS.probes.1.5Kb <- TSS.probes[which(TSS.probes$distance <= 1500),]
######################################################################
mergeMeth.TSS1.5 <- mergeMeth[rownames(mergeMeth)%in%rownames(TSS.probes.1.5Kb),]
mergeCommon.TSS1.5 <- merge(mergeMeth.TSS1.5, geneExp, by="gene")
rownames(mergeCommon.TSS1.5) <- mergeCommon.TSS1.5$Row.names
mergeCommon.TSS1.5$Row.names <- NULL
mergeCommon.TSS1.5[,47:54] <- NULL
mergeCommon.meth.TSS1.5 <- mergeCommon.TSS1.5[,2:46]
#tumorID <- grep("01A", colnames(mergeCommon.meth))
mergeCommon.exp.TSS1.5 <- mergeCommon.TSS1.5[,47:91]
tumorMeth.TSS1.5 <- mergeCommon.meth.TSS1.5[,grep("01A", colnames(mergeCommon.meth.TSS1.5))]
tumorExp.TSS1.5 <- mergeCommon.exp.TSS1.5[,grep("01A", colnames(mergeCommon.exp.TSS1.5))]
colnames(tumorExp.TSS1.5) <- substr(colnames(tumorExp.TSS1.5),1,16)
colnames(tumorMeth.TSS1.5) <- substr(colnames(tumorMeth.TSS1.5),1,16)
#mergeCommon.meth.1.5kb <- mergeCommon[rownames(tumorMeth)%in%rownames(TSS.probes.1.5Kb)]
#tumorMeth.TSS1.5 <- tumorMeth[rownames(tumorMeth)%in%rownames(TSS.probes.1.5Kb),]
#tumorExp.TSS1.5 <- tumorExp[rownames(tumorExp)%in%rownames(TSS.probes.1.5Kb),]

ll <- mapply(function(x,y)cor.test(as.numeric(tumorMeth.TSS1.5[x,]),as.numeric(tumorExp.TSS1.5[y,]), method = "spearman", alternative = "t"),
             1:nrow(tumorMeth.TSS1.5),
             1:nrow(tumorExp.TSS1.5),
             SIMPLIFY=FALSE)
cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
spearman.p <- cbind(cor.value, p.value)
rownames(spearman.p) <- rownames(mergeCommon.TSS1.5)
spearman.p <- data.frame(spearman.p)
spearman.p$gene <- mergeCommon.TSS1.5$gene
#spearman.p$adjP <- p.adjust(spearman.p$p.value, method = c("BH"))
spearman.p1 <- spearman.p[which((spearman.p$p.value <= 0.05)&abs(spearman.p$cor.value) >=0.1),]
spearman.p1 <- spearman.p1[order(spearman.p1$gene),]
spearman.p1.neg <- spearman.p1[which(spearman.p1$cor.value < 0),]
########################
save.image("SpearmanCorrCancer_getNearestTSS.RData")