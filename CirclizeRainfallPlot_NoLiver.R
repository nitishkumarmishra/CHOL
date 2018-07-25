library(gtrellis);library(circlize); library(ComplexHeatmap); library(naturalsort)
setwd("C:/Users/nitish.mishra/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffMeth")
load("DiffMeth.RData")
rm(list= ls()[!(ls() %in% c('results'))])

Liver <- read.csv(file = "Liver_Genes.txt", header = TRUE, sep = "\t")
Liver_ProteinAtlas <- read.csv(file = "ProteinAtlas_tissue_specificity_rna_liver_elevated.tab", header = TRUE, sep = "\t")
Liver_ProteinAtlas_Uniq <- Liver_ProteinAtlas[Liver_ProteinAtlas$RNA.tissue.category!='Group enriched',]
Liver_Gene <- union(Liver$Genes, Liver_ProteinAtlas_Uniq$Gene)

results_Liver <- results[!results$gene%in%Liver$Genes,]


bed <- results_Liver[,c("CHR", "MAPINFO", "MAPINFO", "log2FC")]
bed$CHR <- paste("chr", bed$CHR, sep = "")
### change colname, add start, end, foldchange in dataframe
names(bed)[1] <- c("chr"); names(bed)[2] <- c("start"); names(bed)[3] <- c("end")
bed$direction <- ifelse(bed$log2FC >0, "Hyper", "Hypo")
bed$log2FC <- NULL; 
bed <- bed[naturalorder(bed$chr),]### Order dataframe with chromosome name
DMR_hyper = bed[bed$direction == "Hyper", ]
DMR_hypo = bed[bed$direction == "Hypo", ]
bed_list = list(DMR_hyper, DMR_hypo)
pdf("CirclizeRainfallPlot_NoLiver_logFC2.pdf", width = 8, height = 8)
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
circos.genomicRainfall(bed_list, pch = 16, cex = 0.3, col = c("red", "blue4"))
circos.genomicDensity(DMR_hyper, col = c("red"), track.height = 0.1, window.size = 10000000)
circos.genomicDensity(DMR_hypo, col = c("blue4"), track.height = 0.1, window.size = 10000000)
dev.off()



#### Density plot, density plot by using all dm-CpGs ####
pdf("CirclizeRainfallPlotAllMethDensity_NoLiver_logFC2.pdf", width = 8, height = 8)
circos.initializeWithIdeogram(chromosome.index = paste0("chr", 1:22))
circos.genomicRainfall(bed_list, pch = 16, cex = 0.3, col = c("red", "blue4"))
circos.genomicDensity(bed, col = c("black"), track.height =0.2, window.size = 10000000)
dev.off()
####################################
## Highlight the region with Slice
draw.sector(220, 230, lwd = 2, rou1 = 1, rou2 = NULL, border = "red", clock.wise = FALSE)
################

distance_liver <- rainfallTransform(bed)
density_liver <- genomicDensity(bed, window.size = 1000000)

############# Save differential methylation results #################
results_Liver_0.2 <- results_Liver[abs(results_Liver$MeanDiff) >=0.2,]
results_Liver_0.3 <- results_Liver[abs(results_Liver$MeanDiff) >=0.3,]
write.csv(results_Liver, file = "DiffMethDeltaNoLiver.0.10.logFC2.txt")
write.csv(results_Liver_0.2, file = "DiffMethDeltaNoLiver.0.2.logFC2.txt")
write.csv(results_Liver_0.3, file = "DiffMethDeltaNoLiver.0.3.logFC2.txt")
write.csv(results_Liver_0.2, file = "DiffMethDeltaNoLiver.0.2.logFC2.csv")
write.csv(results_Liver_0.3, file = "DiffMethDeltaNoLiver.0.3.logFC2.csv")

save.image("CirclizeRainfall_NoLiver_logFC2.RData")
