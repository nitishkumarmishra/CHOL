library(DMRcate)
library(missMethyl)
library(minfi)
library(impute)
#getTCGA("CHOL", Meth = TRUE, RNA = TRUE, Clinic = TRUE,   basedir="~", RNAtype = "gene", Methfilter = 0.2)
#load("CHOL_clinic.rda")
load("CHOL_meth.rda")
#load("CHOL_RNA.rda")
####################################
Meth <- rmSNPandCH(Meth, dist = 2, mafcut = 0.05, rmXY = TRUE, rmcrosshyb = TRUE)
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
c1 <- ncol(Tumor.impue.BMIQ)
c2 <- ncol(Normal.impue.BMIQ)
# ####################################
# colnames(Meth) <- substr(colnames(Meth),1,16)
# numNAs <- rowSums(is.na(Meth))
# Meth <- Meth[!(numNAs > dim(Meth)[2]*0.25),]
# Meth.impue <- impute.knn(Meth, k = 15, rowmax = 0.2, colmax = 0.2, rng.seed = 12345)
# Meth.impue <- Meth.impue$data
# ####################################
# myNorm <- champ.norm(beta = Meth.impue, fromIDAT = F, methValue = "B", norm = "BMIQ", fromFile = FALSE, betaFile = Meth.impue, filterXY = FALSE, plotBMIQ = TRUE, resultsDir = paste(getwd(), "resultsChamp", sep = "/"))
# Meth.impue.BMIQ <- myNorm$beta
# Tumor.BMIQ <- Meth.impue.BMIQ[,which(substr(colnames(Meth.impue.BMIQ),14,15)=="01")]
# Normal.BMIQ <- Meth.impue.BMIQ[,which(substr(colnames(Meth.impue.BMIQ),14,15)=="11")]
# BMIQ.Meth <- cbind(Tumor.BMIQ, Normal.BMIQ)
# #BMIQ.Meth.NoSNP <- rmSNPandCH(BMIQ.Meth, dist = 2, mafcut = 0.05)
# BMIQ.Meth <- BMIQ.Meth[-grep("rs", rownames(BMIQ.Meth)),]
# BMIQ.Meth.M <- logit2(BMIQ.Meth)
# c1 <- ncol(Tumor.BMIQ)
# c2 <- ncol(Normal.BMIQ)
###########DMRcate ##################
#BMIQ.Meth.M.NoSNP <- rmSNPandCH(BMIQ.Meth.M, dist = 2, mafcut = 0.05)
#BMIQ.Meth.M.NoSNP.1 <- BMIQ.Meth.M.NoSNP[-grep("rs", rownames(BMIQ.Meth.M.NoSNP)),]
#### Design contrast matrix #########
#design = cbind(Tumor = c(rep(1,  c1), rep(0, c2)), Normal = c(rep(0, c1), rep(1, c2)))
#contrast.matrix = makeContrasts(TumorVsNormal = Tumor-Normal, levels = design)
groups <- as.factor(c(rep("Tumor",c1),rep("Normal",c2)))
design<-model.matrix(~0+groups)
colnames(design)=levels(groups)
contrast.matrix <- makeContrasts(Tumor-Normal, levels = design)
fitvar.contr <- varFit(BMIQ.Meth.M, design=design, coef=c(1,2))
fitvar.contr <- contrasts.varFit(fitvar.contr,contrasts=contrast.matrix)
results <- topVar(fitvar.contr,coef=1, sort = TRUE, number = 1000)
write.csv(results, file = "missMethylDiffVar.csv")
save.image("missMethylDiffVar.RData")
##################
### Pathway analysis by using missMethyl
library(missMethyl)
sigCpGs <- rownames(results.dms)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
gost <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(methFit), collection="GO")
gst <- gometh(sig.cpg=sigCpGs, all.cpg=rownames(methFit), collection="KEGG")
