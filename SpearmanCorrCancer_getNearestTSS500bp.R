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
mergeMeth.TSS500bp <- mergeMeth[rownames(mergeMeth)%in%rownames(TSS.probes.1.5Kb),]
mergeCommon.TSS500bp <- merge(mergeMeth.TSS500bp, geneExp, by="gene")
rownames(mergeCommon.TSS500bp) <- mergeCommon.TSS500bp$Row.names
mergeCommon.TSS500bp$Row.names <- NULL
mergeCommon.TSS500bp[,47:54] <- NULL
mergeCommon.meth.TSS500bp <- mergeCommon.TSS500bp[,2:46]
#tumorID <- grep("01A", colnames(mergeCommon.meth))
mergeCommon.exp.TSS500bp <- mergeCommon.TSS500bp[,47:91]
tumorMeth.TSS500bp <- mergeCommon.meth.TSS500bp[,grep("01A", colnames(mergeCommon.meth.TSS500bp))]
tumorExp.TSS500bp <- mergeCommon.exp.TSS500bp[,grep("01A", colnames(mergeCommon.exp.TSS500bp))]
colnames(tumorExp.TSS500bp) <- substr(colnames(tumorExp.TSS500bp),1,16)
colnames(tumorMeth.TSS500bp) <- substr(colnames(tumorMeth.TSS500bp),1,16)
#mergeCommon.meth.1.5kb <- mergeCommon[rownames(tumorMeth)%in%rownames(TSS.probes.1.5Kb)]
#tumorMeth.TSS500bp <- tumorMeth[rownames(tumorMeth)%in%rownames(TSS.probes.1.5Kb),]
#tumorExp.TSS500bp <- tumorExp[rownames(tumorExp)%in%rownames(TSS.probes.1.5Kb),]

ll <- mapply(function(x,y)cor.test(as.numeric(tumorMeth.TSS500bp[x,]),as.numeric(tumorExp.TSS500bp[y,]), method = "spearman", alternative = "t"),
             1:nrow(tumorMeth.TSS500bp),
             1:nrow(tumorExp.TSS500bp),
             SIMPLIFY=FALSE)
cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
spearman.p <- cbind(cor.value, p.value)
rownames(spearman.p) <- rownames(mergeCommon.TSS500bp)
spearman.p <- data.frame(spearman.p)
spearman.p$gene <- mergeCommon.TSS500bp$gene
#spearman.p$adjP <- p.adjust(spearman.p$p.value, method = c("BH"))
spearman.p1 <- spearman.p[which((spearman.p$p.value <= 0.05)&abs(spearman.p$cor.value) >=0.2),]
spearman.p1 <- spearman.p1[order(spearman.p1$gene),]
spearman.p1.neg <- spearman.p1[which(spearman.p1$cor.value < 0),]
########################
save.image("SpearmanCorrCancer_getNearestTSS500bp.RData")
########################################
#Pearson's correlation nalysis
ll <- mapply(function(x,y)cor.test(as.numeric(tumorMeth.TSS500bp[x,]),as.numeric(tumorExp.TSS500bp[y,]), method = "pearson", alternative = "t"),
             1:nrow(tumorMeth.TSS500bp),
             1:nrow(tumorExp.TSS500bp),
             SIMPLIFY=FALSE)
cor.value <- sapply(ll,'[[','estimate')
p.value <- sapply(ll,'[[','p.value')
pearson.p <- cbind(cor.value, p.value)
rownames(pearson.p) <- rownames(mergeCommon.TSS500bp)
pearson.p <- data.frame(pearson.p)
pearson.p$gene <- mergeCommon.TSS500bp$gene
#pearson.p$adjP <- p.adjust(pearson.p$p.value, method = c("BH"))
pearson.p1 <- pearson.p[which((pearson.p$p.value <= 0.05)&abs(pearson.p$cor.value) >=0.2),]
pearson.p1 <- pearson.p1[order(pearson.p1$gene),]
pearson.p1.neg <- pearson.p1[which(pearson.p1$cor.value < 0),]
########################
save.image("PearsonCorrCancer_getNearestTSS500bp.RData")