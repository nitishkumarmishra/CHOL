library(ggpubr)
load("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffMeth/Survival_analysis_Promoter_Methylation_Hyper.RData")

tmp <- as.data.frame(t(BMIQ.Meth.Hyper.TSS))
rownames(tmp) <- substr(rownames(tmp), 1, 16)
tmp.select <- tmp %>% select("cg27362525", "cg26597242")
#colnames(tmp.select) <- c("Mir551b", "Mir22")
#tmp.select$Mir551b <- log2(tmp.select$Mir551b+1)
#tmp.select$Mir22 <- log(tmp.select$Mir22+1)/2
tmp.select$Sample <- ifelse(substr(rownames(tmp.select), 14,15)=="01", "Tumor", "Normal")



p <- ggboxplot(tmp.select, x = "Sample", y = "cg27362525", color = "Sample",
               palette = c( "red4","blue4"), add = "jitter", ylab = "DNA methylation Beta value for cg27362525")
p + stat_compare_means(method = "t.test", label.y = max(tmp.select$cg27362525)+0.5)
dev.print(pdf, 'cg27362525 Boxplot.pdf ', width = 6, height = 6)



p <- ggboxplot(tmp.select, x = "Sample", y = "cg26597242", color = "Sample",
               palette = c( "red4","blue4"), add = "jitter", ylab = "DNA methylation Beta value for cg26597242")
p + stat_compare_means(method = "t.test", label.y = max(tmp.select$cg26597242)+0.5)
dev.print(pdf, 'cg26597242 Boxplot.pdf ', width = 6, height = 6)

