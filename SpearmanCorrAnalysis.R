### R code to correlate gene expression and DNA methylation

#DiffMeth <- read.csv("DiffMethDelta.0.10.txt", sep = "\t")## Probe with Delta beta 0.2
#DiffExp <- read.csv("mRNA_results.0.01.txt", row.names = 1)## Diff exp genes with adjP < 0.01
#DiffExp$gene <- rownames(DiffExp)
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
##########Correlation for only TSS200/TSS1500
mergeMeth.TSS <- mergeMeth[grep("TSS",mergeMeth$feature),]
#mergeMeth.1 <- mergeMeth[DiffMeth$X,]
#mergeExp.1 <- mergeExp[rownames(DiffExp),]
#mergeMeth.1$Row.names <- NULL
mergeCommon <- merge(mergeMeth.TSS, geneExp, by="gene")
rownames(mergeCommon) <- mergeCommon$Row.names
mergeCommon$Row.names <- NULL
mergeCommon[,47:54] <- NULL
mergeCommon.meth <- mergeCommon[,2:46]
#tumorID <- grep("01A", colnames(mergeCommon.meth))
mergeCommon.exp <- mergeCommon[,47:91]
tumorMeth <- mergeCommon.meth[,grep("01A", colnames(mergeCommon.meth))]
tumorExp <- mergeCommon.exp[,grep("01A", colnames(mergeCommon.exp))]
colnames(tumorExp) <- substr(colnames(tumorExp),1,16)
colnames(tumorMeth) <- substr(colnames(tumorMeth),1,16)
################# Spearman's correlation analysis ########
ll <- mapply(function(x,y)cor.test(as.numeric(mergeCommon.meth[x,]),as.numeric(mergeCommon.exp[y,]), method = "spearman", alternative = "t"),
             1:nrow(mergeCommon.meth),
             1:nrow(mergeCommon.exp),
             SIMPLIFY=FALSE)
cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
spearman.p <- cbind(cor.value, p.value)
rownames(spearman.p) <- rownames(mergeCommon.meth)
spearman.p <- data.frame(spearman.p)
spearman.p$gene <- mergeCommon$gene
#spearman.p$adjP <- p.adjust(spearman.p$p.value, method = c("BH"))
spearman.p1 <- spearman.p[which((spearman.p$p.value <= 0.05)&abs(spearman.p$cor.value) >=0.1),]
spearman.p1 <- spearman.p1[order(spearman.p1$gene),]
spearman.p1.neg <- spearman.p1[which(spearman.p1$cor.value < 0),]
###################### Save image ########################
save.image(file = "SpearmanCorrTSS.RData")
##########################################################
ll <- mapply(function(x,y)cor.test(as.numeric(tumorMeth[x,]),as.numeric(tumorExp[y,]), method = "spearman", alternative = "t"),
             1:nrow(tumorMeth),
             1:nrow(tumorExp),
             SIMPLIFY=FALSE)
cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
spearman.p <- cbind(cor.value, p.value)
rownames(spearman.p) <- rownames(mergeCommon.meth)
spearman.p <- data.frame(spearman.p)
spearman.p$gene <- mergeCommon$gene
#spearman.p$adjP <- p.adjust(spearman.p$p.value, method = c("BH"))
spearman.p1 <- spearman.p[which((spearman.p$p.value <= 0.05)&abs(spearman.p$cor.value) >=0.1),]
spearman.p1 <- spearman.p1[order(spearman.p1$gene),]
spearman.p1.neg <- spearman.p1[which(spearman.p1$cor.value < 0),]
###################### Save image ########################
save.image(file = "SpearmanCorrCancerTSS.RData")
##########################################################

