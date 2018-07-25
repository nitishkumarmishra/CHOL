library(pacman)
### p_load command can load several package at a time
p_load(ELMER, TCGAbiolinks, impute, matrixStats, ChAMP, limma, multtest, DESeq2, edgeR, DMRcate, missMethyl)
data(dmrcatedata)
#p_loaded() ### To check the loaded packages
#p_unload(negate = TRUE)## Unload all packages
############# load TCGA level3 data ############
#getTCGA("CHOL", Meth = TRUE, RNA = TRUE, Clinic = TRUE,   basedir="~", RNAtype = "gene", Methfilter = 0.2)
#load("CHOL_clinic.rda")
load("CHOL_meth.rda")
#load("CHOL_RNA.rda")
colnames(Meth) <- substr(colnames(Meth), 1, 16)
####################################
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
BMIQ.Meth <- rmSNPandCH(BMIQ.Meth, dist = 1, mafcut = 0.05, rmXY = TRUE, rmcrosshyb = TRUE)
#BMIQ.Meth <- BMIQ.Meth[-grep("rs", rownames(BMIQ.Meth)),]## Remove rs probes
colnames(BMIQ.Meth) <- substr(colnames(BMIQ.Meth), 6, 16)
BMIQ.Meth.M <- logit2(BMIQ.Meth)
c1 <- ncol(Tumor.impue.BMIQ)
c2 <- ncol(Normal.impue.BMIQ)
###########DMRcate ##################
#BMIQ.Meth.M.NoSNP <- rmSNPandCH(BMIQ.Meth.M, dist = 2, mafcut = 0.05, rmXY = TRUE, rmcrosshyb = TRUE)
#BMIQ.Meth.M.NoSNP.1 <- BMIQ.Meth.M.NoSNP[-grep("rs", rownames(BMIQ.Meth.M.NoSNP)),]
#### Design contrast matrix #########
#design = cbind(Tumor = c(rep(1,  c1), rep(0, c2)), Normal = c(rep(0, c1), rep(1, c2)))
#contrast.matrix = makeContrasts(TumorVsNormal = Tumor-Normal, levels = design)
groups <- as.factor(c(rep("Tumor",c1),rep("Normal",c2)))
design<-model.matrix(~0+groups)
colnames(design)=levels(groups)
contrast.matrix <- makeContrasts(Tumor-Normal, levels = design)
######### DMRcate anaysis ###########
#myannotation <- cpg.annotate("array",BMIQ.Meth.M, analysis.type = "differential", design = design, contrasts = TRUE, cont.matrix = contrast.matrix, coef = "Tumor - Normal", fdr = 0.01)
myannotation <- cpg.annotate("array",BMIQ.Meth.M, analysis.type = "differential", design = design, contrasts = TRUE, cont.matrix = contrast.matrix, coef = "Tumor - Normal")
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2, p.adjust.method = "BH")
#### Annotate overlapping promoter regions (+/- 2000 bp from TSS)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
groups <- c(Tumor="magenta", Normal="forestgreen")
type <- as.factor(c(rep("Tumor", c1), rep("Normal", c2)))
cols <- groups[as.character(type)]
samps <- c(1:36, 36+(1:9))
#DMR.plot(ranges=results.ranges, dmr=1, CpGs=BMIQ.Meth, phen.col=cols, genome="hg19", samps=samps)
par(mfrow=c(1,1))
#DMR.plot(ranges=results.ranges, dmr=1, CpGs=BMIQ.Meth, phen.col=cols, genome="hg19", samps=samps, showSampleNames = TRUE, cex.sampleNames = 0.8, separator = 1, pch=18, toscale=TRUE, plotmedians=TRUE,legend = TRUE)
DMR.plot(ranges=results.ranges, dmr=1, CpGs=BMIQ.Meth, phen.col=cols, genome="hg19", samps=samps, showSampleNames = FALSE, col.axis="Black", col.title="Black",cex.sampleNames = 0.8, col.sampleNames="Black", 
         separator = 1.75, pch=18, toscale=TRUE, plotmedians=TRUE,legend = TRUE, cex = 0.7, fontcolor="Black", fontcolor.group="Black", fontsize.group=12, fontcolor.feature = 1, cex.feature = 0.7, min.width = 3, min.distance = 5, cex.title=0.6, rotation.title=360)
rm(c1, c2, numNAs, myNorm, Meth, tx.hg19, tx.hg38, tx.mm10, probe.features, probeInfoALL.lv)

############## Conver Grange file in matrix #######
#### List of DMR by using DMRcate #################
df1 <- data.frame(seqnames=seqnames(results.ranges),
                  starts=start(results.ranges),
                  ends=end(results.ranges),
                  names=c(rep(".", length(results.ranges))),
                  strand=strand(results.ranges)
)
df <- mcols(results.ranges)
df2 <- cbind(df1, df)
write.csv(df2, file = "DMRcateDMR.csv")
save.image("DMRcate.RData")
