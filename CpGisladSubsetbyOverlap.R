### DMRcate result and other data
coords <- strsplit2(dmrcoutput$results$coord, ":")
chrom <- coords[,1]
start <- sapply(strsplit(as.character(coords[,2]), "-"), "[[", 1)
end <- sapply(strsplit(as.character(coords[,2]), "-"), "[[", 2)
start <- as.numeric(start); end <- as.numeric(end)
coordsCpG <- GRanges(seqnames=Rle(chrom),
                      ranges=IRanges(start=start, end=end),
                      strand=Rle(strand(rep("*",nrow(coords)))))

islandHMM = read.csv(paste("../model-based-cpg-islands-hg19.txt", sep="/"), sep="\t", stringsAsFactors=FALSE, header=TRUE)
islandData <- GRanges(seqnames=Rle(islandHMM$chr),
                      ranges=IRanges(start=islandHMM$start, end=islandHMM$end),
                      strand=Rle(strand(rep("*",nrow(islandHMM)))))

dnase <- read.csv(paste("../wgEncodeRegDnaseClusteredV3.bed",
                        sep="/"), sep="\t",stringsAsFactors=FALSE,header=FALSE)
dnaseData <- GRanges(seqnames=dnase[,1],
                     ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                     strand=Rle(rep("*",nrow(dnase))), data=dnase[,5])
### Mean CpG length=608.35,, Avg No. of CpGs in per island = 48.44
# Mean DMR length = 1012.21, average number of Cps in per DMR=6.92
## Average 146 base pair for 1 CpGs in DMR
#DMRinCpG <- subsetByOverlaps(results.ranges, islandData, minoverlap = 152)## 25% of avg CpG length
### DMR in CpG islands
DMRinCpG <- subsetByOverlaps(results.ranges, islandData)## DMR in CpG
DMRinCpGandDNAse <- subsetByOverlaps(DMRinCpG, dnaseData)## DMR in DNAse region
df1 <- data.frame(seqnames=seqnames(DMRinCpG),
           starts=start(DMRinCpG),
           ends=end(DMRinCpG),
           names=c(rep(".", length(DMRinCpG))),
           strand=strand(DMRinCpG)
)
df <- mcols(DMRinCpG)
df2 <- cbind(df1, df)
write.csv(df2, file = "CpGsDMRcate.csv")
### DMR in CpG islands e hypersensitive region
df1 <- data.frame(seqnames=seqnames(DMRinCpGandDNAse),
                  starts=start(DMRinCpGandDNAse),
                  ends=end(DMRinCpGandDNAse),
                  names=c(rep(".", length(DMRinCpGandDNAse))),
                  strand=strand(DMRinCpGandDNAse)
)
df <- mcols(DMRinCpGandDNAse)
df2 <- cbind(df1, df)
write.csv(df2, file = "DMRinCpGandDNAse.csv")




save.image("CpGsDMRcate.RData")
