## R code to get probes with 1.5 kb from TSS
## By using getNearestTSS from FDb.InfiniumMethylation.hg19
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
save.image(file = "Nitish_TSS_1500bp.RData")
