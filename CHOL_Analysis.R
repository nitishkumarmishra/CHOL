library(pacman)
### p_load command can load several package at a time
p_load(ChAMP, DMRcate, ELMER, impute, TCGAbiolinks, impute, matrixStats, limma, multtest, DESeq2, edgeR)
#p_loaded() ### To check the loaded packages
#p_unload(negate = TRUE)## Unload all packages
############# load TCGA level3 data ############
getTCGA("CHOL", Meth = TRUE, RNA = TRUE, Clinic = TRUE,   basedir="~", RNAtype = "gene", Methfilter = 0.2)
load("CHOL_clinic.rda")
load("CHOL_meth.rda")
load("CHOL_RNA.rda")
colnames(Meth) <- substr(colnames(Meth),1,16)
####################################
####################################
Meth <- rmSNPandCH(Meth, dist = 10, mafcut = 0.05, rmXY = TRUE, rmcrosshyb = TRUE)
Meth <- Meth[-grep("rs", rownames(Meth)),]## Remove rs probes
Tumor.Meth <- Meth[,which(substr(colnames(Meth),14,15)=="01")]
Normal.Meth <- Meth[,which(substr(colnames(Meth),14,15)=="11")]
numNAs <- rowSums(is.na(Tumor.Meth))
Tumor.Meth <- Tumor.Meth[!(numNAs > dim(Tumor.Meth)[2]*0.25),]## Remove probe which have >25% NAs
numNAs <- rowSums(is.na(Normal.Meth))
Normal.Meth <- Normal.Meth[!(numNAs > dim(Normal.Meth)[2]*0.25),]
commonID <- intersect(rownames(Tumor.Meth), rownames(Normal.Meth))
Tumor.Meth <- Tumor.Meth[commonID,]
Normal.Meth <- Normal.Meth[commonID,]
Tumor.impue <- impute.knn(Tumor.Meth, k = 15, rowmax = 0.25, colmax = 0.25, rng.seed = 12345)
Tumor.impue <- Tumor.impue$data
Normal.impue <- impute.knn(Normal.Meth, k = 15, rowmax = 0.25, colmax = 0.25, rng.seed = 12345)
Normal.impue <- Normal.impue$data
myNorm <- champ.norm(beta = Tumor.impue, fromIDAT = F, methValue = "B", norm = "BMIQ", fromFile = FALSE, betaFile = Tumor.impue, filterXY = TRUE, plotBMIQ = TRUE, resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
Tumor.impue.BMIQ <- myNorm$beta
myNorm <- champ.norm(beta = Normal.impue, fromIDAT = F, methValue = "B", norm = "BMIQ", fromFile = FALSE, betaFile = Normal.impue, filterXY = TRUE, plotBMIQ = TRUE, resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
Normal.impue.BMIQ <- myNorm$beta
BMIQ.Meth <- cbind(Tumor.impue.BMIQ, Normal.impue.BMIQ)
BMIQ.Meth.M <- logit2(BMIQ.Meth)
#Probeinfo.merge <- merge(Probeinfo, probe.features, by = "row.names") ## Merg
####################################
####################################
# numNAs <- rowSums(is.na(Meth))
# Meth <- Meth[!(numNAs > dim(Meth)[2]*0.25),]
# Meth.impue <- impute.knn(Meth, k = 15, rowmax = 0.2, colmax = 0.2, rng.seed = 12345)
# Meth.impue <- Meth.impue$data
# ####################################
# myNorm <- champ.norm(beta = Meth.impue, fromIDAT = F, methValue = "B", norm = "BMIQ", fromFile = FALSE, betaFile = Meth.impue, filterXY = FALSE, plotBMIQ = TRUE, resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
# Meth.impue.BMIQ <- myNorm$beta
# Tumor.BMIQ <- Meth.impue.BMIQ[,which(substr(colnames(Meth.impue.BMIQ),14,15)=="01")]
# Normal.BMIQ <- Meth.impue.BMIQ[,which(substr(colnames(Meth.impue.BMIQ),14,15)=="11")]
#BMIQ.Meth <- cbind(Tumor.BMIQ, Normal.BMIQ)
#BMIQ.Meth.M <- logit2(BMIQ.Meth)
#####################################
#c1 <- ncol(Tumor.BMIQ)
#c2 <- ncol(Normal.BMIQ)
c1 <- ncol(Tumor.impue.BMIQ)
c2 <- ncol(Normal.impue.BMIQ)
design = cbind(Tumor = c(rep(1,  c1), rep(0, c2)), Normal = c(rep(0, c1), rep(1, c2)))
methFit = lmFit(BMIQ.Meth.M, design)
contrast.matrix = makeContrasts(TumorVsNormal = Tumor-Normal, levels = design)
methFit2 = contrasts.fit(methFit, contrast.matrix)
methFit2 = eBayes(methFit2)
methRresults = topTable(methFit2, coef = 1, number = Inf, p.value = 0.01, sort.by = "p", resort.by="logFC", adjust.method = "BH")
#####################################
## Annotation of limma result file ##
#t.tumor <- as.data.frame(cbind(rowMeans(Tumor.BMIQ), rowMeans(Normal.BMIQ), rowMeans(Tumor.BMIQ) - rowMeans(Normal.BMIQ), rowMeans(Tumor.BMIQ) / rowMeans(Normal.BMIQ), log2(rowMeans(Tumor.BMIQ) / rowMeans(Normal.BMIQ))))
t.tumor <- as.data.frame(cbind(rowMeans(Tumor.impue.BMIQ), rowMeans(Normal.impue.BMIQ), rowMeans(Tumor.impue.BMIQ) - rowMeans(Normal.impue.BMIQ), rowMeans(Tumor.impue.BMIQ) / rowMeans(Normal.impue.BMIQ), log2(rowMeans(Tumor.impue.BMIQ) / rowMeans(Normal.impue.BMIQ))))
colnames(t.tumor) <- c("Tumor", "Normal", "MeanDiff", "FoldChange", "log2FC")
limma.results <- merge(methRresults, t.tumor, by="row.names")
rownames(limma.results) <- limma.results$Row.names
limma.results$Row.names <- NULL
results <- merge(limma.results, probe.features, by="row.names")
results <- results[which(results$adj.P.Val <= 0.01 & abs(results$MeanDiff) >= 0.10 & abs(results$logFC) > 2),]
#results <- results[which(results$adj.P.Val <= 0.01 & abs(results$MeanDiff) >= 0.10 & (results$FoldChange >= 2|results$FoldChange <= 0.5)),]
rownames(results) <- results$Row.names
results$Row.names<- NULL
write.csv(results, file = "DiffMethDelta.0.10.logFC2.csv", sep = "\t")
save.image("DiffMeth_logFC2.RData")
#####################################
### Differential expression #########

