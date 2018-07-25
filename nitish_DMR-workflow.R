library(DMRcate)
a <- read.csv("C:/Users/pet22q/Downloads/BMIQ.CSV")
rownames(a) <- a$X
a <- a[,-1]
betas <- data.matrix(a)

BMIQ.Meth.M <- logit2(betas)
c1 <- length(grep("Tumor", colnames(BMIQ.Meth.M)))
c2 <- length(grep("Normal", colnames(BMIQ.Meth.M)))

groups <- as.factor(c(rep("Tumor",c1),rep("Normal",c2)))
design<-model.matrix(~0+groups)
colnames(design)=levels(groups)

library(limma)
contrast.matrix <- makeContrasts(Tumor - Normal, levels = design)

myannotation <- cpg.annotate("array",BMIQ.Meth.M, analysis.type = "differential", design = design, contrasts = TRUE, cont.matrix = contrast.matrix, coef = "Tumor - Normal")
data(dmrcatedata)
dmrcoutput <- dmrcate(myannotation)
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
groups <- c(Tumor="magenta", Normal="forestgreen")
type <- c(rep("Tumor",c1), rep("Normal", c2))
cols <- groups[as.character(type)]
samps <- c(1:36, 36+(1:9))
DMR.plot(ranges=results.ranges, dmr=1, CpGs=betas, phen.col=cols, genome="hg19", samps=samps)
