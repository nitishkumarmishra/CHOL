### DiffMeth is differentially methylated file
library(gtrellis);library(circlize); library(ComplexHeatmap)
bed <- DiffMeth[,c("CHR", "MAPINFO",  "MAPINFO", "log2FC")]
bed$CHR <- paste("chr", bed$CHR, sep = "")
colnames(bed)[1] <- c("chr")
colnames(bed)[2] <- c("start")
colnames(bed)[3] <- c("end")

pdf("DiffMethRainfallPlot-Delta10.pdf", width = 10, height = 10)
gtrellis_layout(ncol = 4, byrow = FALSE,
                track_axis = TRUE, 
                track_ylim = range(bed[[4]]),
                track_ylab = c("log2FC"),
                add_name_track = TRUE, add_ideogram_track = TRUE)
add_points_track(bed, bed[[4]], pch = 16, size = unit(1, "mm"), gp = gpar(col = ifelse(bed[[4]] > 0, "red", "green")))
dev.off()


####################
library(quantsmooth)
bed$chr <- gsub("chr", "", bed$chr)## remove "chr" from CHR name
annot<-bed[,c("chr", "start")] # chromosome and position
colnames(annot)<-c("CHR","MapInfo")
annot$CHR[annot$CHR==23]<-"X"
annot$CHR[annot$CHR==24]<-"Y"
# annot <- annot[,-(annot$CHR==23)]
# annot <- annot[,-(annot$CHR==24)]
bmp("ChromosomeDistribution.jpg", height=1000,width=1000)
out<-prepareGenomePlot(annot,sexChromosomes=F,organism="hsa",units="bases",paintCytobands=T)
out<-data.frame(out)
logp<- bed$log2FC
x<-out$MapInfo[bed$log2FC<0]
y<-(out$CHR+(logp/10))[bed$log2FC<0]
points(x=x,y=y,pch=20,cex=0.8,col="blue4")
x<-out$MapInfo[bed$log2FC>0]
y<-(out$CHR+(logp/10))[bed$log2FC>0]
points(x=x,y=y,pch=20,cex=0.8,col="red4")
dev.off()
