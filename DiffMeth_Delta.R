library(pacman)
### p_load command can load several package at a time
#p_load(ChAMP, DMRcate, ELMER, impute, TCGAbiolinks, impute, matrixStats, limma, multtest, DESeq2, edgeR)
p_load(ChAMP, DMRcate, impute, matrixStats, limma, multtest, lumi)
#p_loaded() ### To check the loaded packages
#p_unload(negate = TRUE)## Unload all packages
############# load TCGA level3 data ############
#getTCGA("CHOL", Meth = TRUE, RNA = TRUE, Clinic = TRUE,   basedir="~", RNAtype = "gene", Methfilter = 0.2)
#load("CHOL_clinic.rda")
#load("CHOL_meth.rda")
#load("CHOL_RNA.rda")
load("probe.features.Rda")
load("probeInfo_feature_distal.rda")
Meth <- read.csv(file = "CHOL.methylation__humanmethylation450__jhu_usc_edu__Level_3__within_bioassay_data_set_function__data.data.txt", header = TRUE, sep = "\t", stringsAsFactors=FALSE, row.names = 1, check.names = FALSE)
Meth <- Meth[,which(Meth[1,]=='Beta_value')]
Meth <- Meth[-1,]
colnames(Meth) <- gsub("\\.", "-", colnames(Meth))
colnames(Meth) <- substr(colnames(Meth),1,16)

save(Meth, file = "CHOL_Firehose_Beta.rda")
####################################
####################################
Meth1 <- Meth ## JUST FOR BACKUP
Meth <- as.matrix(Meth)
numNAs <- rowSums(is.na(Meth))
Meth <- Meth[!(numNAs > dim(Meth)[2]*0.25),]
Meth <- Meth[-grep("rs", rownames(Meth)),]
Meth <- rmSNPandCH(Meth.impue, dist = 10, mafcut = 0.05, rmXY = TRUE, rmcrosshyb = TRUE)

Meth.impue <- impute.knn(Meth, k = 15, rowmax = 0.2, colmax = 0.2, rng.seed = 12345)
Meth.impue <- Meth.impue$data

# ####################################
#myNorm <- champ.norm(beta = Meth.impue, fromIDAT = F, methValue = "B", norm = "BMIQ", fromFile = FALSE, betaFile = Meth.impue, filterXY = FALSE, plotBMIQ = TRUE, resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
## champ.norm command changed in latest champ 
myNorm <- champ.norm(beta = Meth.impue, method = "BMIQ", plotBMIQ = FALSE, arraytype="450K", resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
#Meth.impue.BMIQ <- myNorm$beta
Meth.impue.BMIQ <- myNorm ## This also changed
Tumor.BMIQ <- Meth.impue.BMIQ[,which(substr(colnames(Meth.impue.BMIQ),14,15)=="01")]
Normal.BMIQ <- Meth.impue.BMIQ[,which(substr(colnames(Meth.impue.BMIQ),14,15)=="11")]
BMIQ.Meth <- cbind(Tumor.BMIQ, Normal.BMIQ)
#BMIQ.Meth.M <- logit2(BMIQ.Meth)## minfi logit2
BMIQ.Meth.M <- beta2m(BMIQ.Meth)# lumi beta2m, but both are same
#####################################
c1 <- ncol(Tumor.BMIQ)
c2 <- ncol(Normal.BMIQ)

design <- model.matrix(~0 + factor(c(rep(2, c1), rep(1, c2))))
colnames(design) <- c("Normal", "Tumor")
cont.matrix <- makeContrasts("Tumor-Normal", levels = design)
#design = cbind(Tumor = c(rep(1,  c1), rep(0, c2)), Normal = c(rep(0, c1), rep(1, c2)))
methFit = lmFit(BMIQ.Meth, design)
#cont.matrix = makeContrasts(TumorVsNormal = Tumor-Normal, levels = design)
methFit2 = contrasts.fit(methFit, cont.matrix)
methFit2 = eBayes(methFit2)
#methRresults <- topTable(methFit2, adjust.method = "BH", number = length(methFit2), sort.by = "p", resort.by="logFC")
methRresults = topTable(methFit2, number = Inf, p.value = 0.01, sort.by = "p", resort.by="logFC", adjust.method = "BH")
#####################################
## Annotation of limma result file ##
#t.tumor <- as.data.frame(cbind(rowMeans(Tumor.BMIQ), rowMeans(Normal.BMIQ), rowMeans(Tumor.BMIQ) - rowMeans(Normal.BMIQ), rowMeans(Tumor.BMIQ) / rowMeans(Normal.BMIQ), log2(rowMeans(Tumor.BMIQ) / rowMeans(Normal.BMIQ))))
t.tumor <- as.data.frame(cbind(rowMeans(Tumor.BMIQ), rowMeans(Normal.BMIQ), rowMeans(Tumor.BMIQ) - rowMeans(Normal.BMIQ), rowMeans(Tumor.BMIQ) / rowMeans(Normal.BMIQ), log2(rowMeans(Tumor.BMIQ) / rowMeans(Normal.BMIQ))))
colnames(t.tumor) <- c("Tumor", "Normal", "MeanDiff", "FoldChange", "log2FC")
limma.results <- merge(methRresults, t.tumor, by="row.names")
rownames(limma.results) <- limma.results$Row.names
limma.results$Row.names <- NULL
results <- merge(limma.results, probe.features, by="row.names")
#results.dms <- results[which(results$adj.P.Val <= 0.01 & abs(results$MeanDiff) >= 0.20 & abs(results$FoldChange) > 2),]
results.dms <- results[which(results$adj.P.Val <= 0.005 & abs(results$MeanDiff) >= 0.20),]
results.dms <- results.dms[order(results.dms$MeanDiff, decreasing = TRUE),]
#results <- results[which(results$adj.P.Val <= 0.01 & abs(results$MeanDiff) >= 0.10 & (results$FoldChange >= 2|results$FoldChange <= 0.5)),]
rownames(results.dms) <- results.dms$Row.names
results.dms$Row.names<- NULL
results.dms.hyper <- results.dms[results.dms$logFC >0,]
results.dms.hypo <- results.dms[results.dms$logFC <0,]
write.csv(results.dms, file = "DiffMethDelta.0.20.csv")

#####################################
TSS1500.gene <- read.csv("TSS1500.gene", header = TRUE, sep = "\t", row.names = 1)
length(intersect(rownames(TSS1500.gene), rownames(results.dms.hyper)))



results.dms.01 <- results[which(results$adj.P.Val <= 0.01 & abs(results$MeanDiff) >= 0.20),]
results.dms.hyper.01 <- results.dms.01[results.dms.01$logFC >0,]
results.dms.hypo.01 <- results.dms.01[results.dms.01$logFC <0,]

### Differential expression #########
lincRNA.CpG <- read.csv(file = "lincRNACpG.txt", header = TRUE, sep = "\t")
results.dms.lincRNA <- lincRNA.CpG[lincRNA.CpG$Illumina.ID%in%rownames(results.dms),]
rownames(results.dms.lincRNA) <- results.dms.lincRNA$Illumina.ID
write.csv(results.dms.lincRNA.FC, file = "results.dms.lincRNA.txt")
results.dms.lincRNA.FC <- results.dms[results.dms.lincRNA$Illumina.ID,]
dim(results.dms.lincRNA.FC[results.dms.lincRNA.FC$log2FC < 0 ,])


#################################
library(missMethyl)
sigCpGs <- rownames(results.dms)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
gost <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(methFit), collection="GO")
gost <- gost[order(gost$FDR, decreasing = FALSE),]
gst <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(methFit), collection="KEGG")
gst <- gst[order(gst$FDR, decreasing = FALSE),]
gost <- gost[gost$FDR < 0.01,]
gst <- gst[gst$FDR < 0.01,]
write.csv(gost, file = "DiffMeth_Delta0.2_GO.csv")
write.csv(gst, file = "DiffMeth_Delta0.2_KEGG.csv")
##############################
save.image("DiffMeth_Delta0.2.RData")
