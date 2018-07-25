library(ComplexHeatmap)
library(circlize)
#load("ComplexHeatmap_DEG.RData")
load("C:/Users/nitish.mishra/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffExp/DiffExp_Comparative_analysis.RData")
rm(list= ls()[!(ls() %in% c('cnts', 'cnts.liver', 'common_names', 'list_of_data'))])
cancerID <- grep("01A", colnames(cnts.liver))
normalID <- grep("11A", colnames(cnts.liver))
Meth <- log(cnts.liver+1)
colnames(Meth) <- substr(colnames(Meth), 1, 16)
rownames(Meth) <- sapply(strsplit(rownames(Meth),"\\|"),'[[',1)
gene <- sapply(strsplit(common_names,"\\|"),'[[',1)### Common from Limma_voom, edgeR and DESeq2 all
Meth1 <- Meth <- Meth[gene,]


#####################
load("C:/Users/nitish.mishra/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/CHOL_clinic.rda")
gender <- Clinic$gender
race <- Clinic$race
histology <- Clinic$histologic_diagnosis
histology <- gsub("Cholangiocarcinoma; ", "", histology)
grade <- Clinic$tumor_grade
barcode <- Clinic$bcr_patient_barcode
Status <- race #### I have to put Status on top of plot, So it should be fisrst column id data
data <- data.frame(Status=c(Status), Gender=c(gender), Race=c(race), Histology=c(histology), Grade=c(grade))
rownames(data) <- barcode
Ids <- substr(colnames(Meth1), 1, 12)
data <- data[c(Ids),]
data$Status <- c(rep("Tumor", length(cancerID)), rep("Normal", length(normalID)))
data$Histology <- gsub("hilar/perihilar", "Hillar", data$Histology)
data$Histology <- gsub("intrahepatic", "Intrahepatic", data$Histology)
data$Histology <- gsub("distal", "Distal", data$Histology)
data$Race <- gsub("BLACK OR AFRICAN AMERICAN", "Black", data$Race)
data$Race <- gsub("ASIAN", "Asian", data$Race)
data$Race <- gsub("WHITE", "White", data$Race)
data$Gender <- gsub("^MALE$", "Male", data$Gender)
data$Gender <- gsub("FEMALE", "Female", data$Gender)

results <- list_of_data[[2]]
rownames(results) <- sapply(strsplit(rownames(results),"\\|"),'[[',1)
results <- results[gene,]
foldchange <- results$log2FoldChange

set.seed(12345)
pdf("Heatmap2.pdf", width = 12, height = 12)
df <- data
ha <- HeatmapAnnotation(df = df, col = list(Gender= c("Male" = "red4", "Female"="blue4"),
                                            Status = c("Tumor" = "red4", "Normal"="blue4"), 
                                            Race= c("Asian" = "blue4", "White"="red4", "Black"="green4"),
                                            Histology= c("Intrahepatic" = "red4", "Hillar"="green4", "Distal"="blue4"), 
                                            Grade= c("G1" = "blue4", "G2"="green4", "G3"="red4", "G4"="gray10")), gap = unit(c(5,2,2,2,2,2,2,2), "mm"), show_legend = TRUE)
ht = Heatmap(Meth1, clustering_method_columns = "ward.D2", name = "Expression", col = colorRamp2(c(0, 4, 8, 16), c("green4","mediumseagreen","red", "red4")),
             cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 5, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), row_dend_side  = "right", combined_name_fun = NULL, column_title_gp = gpar(fontsize = 8)) +
  #Heatmap(AUC, name = "AUC", col = colorRamp2(c(0.75, 0.88, 1.0), c("blue4", "purple", "red4")), column_names_gp = gpar(fontsize = 8))+
  Heatmap(foldchange, name = "logFC", show_row_names = FALSE, col = colorRamp2(c(11, 0, -5, -10), c("red4","red" ,"blue1", "blue4")), column_names_gp = gpar(fontsize = 8))
draw(ht,heatmap_legend_side = "left", column_title_gp = gpar(fontsize = 12, fontface = "bold"))
for(an in colnames(df)) {
  decorate_annotation(an, {
    grid.text(an, unit(1, "npc") + unit(2, "mm"), 0.5, default.units = "npc", just = "left")
  })}
dev.off()


save.image(file = "ComplexHeatmap_DEG_New.RData")
################################################
set.seed(1110)
ha = HeatmapAnnotation(df = data.frame(type = c(rep("Chol", length(cancerID)), rep("Normal", length(normalID)))), col = list(type = c("Chol" = "red4", "Normal"="blue4")))
ht_list = Heatmap(Meth1, clustering_method_columns = "ward.D2", name = "Expression", col = colorRamp2(c(0, 8, 16), c("green4","red", "red4")),
                  cluster_columns = TRUE, cluster_rows = TRUE, show_column_names = TRUE, top_annotation = ha, column_names_gp = gpar(fontsize = 8), km = 6, show_row_names = FALSE, row_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 8)) +
  Heatmap(foldchange, name = "Fold Change", show_row_names = FALSE, col = colorRamp2(c(11, 0, -5, -10), c("red4","red" ,"blue1", "blue4")), column_names_gp = gpar(fontsize = 8))
pdf("Chol_Exp_All.pdf", width = 10, height = 10)
draw(ht_list, newpage = FALSE, column_title = "Comprehensive differential gene expression analysis", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
dev.off()