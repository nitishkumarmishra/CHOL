library(pacman)
p_load('TCGAbiolinks', 'survival', 'survminer')

setwd("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffMeth/")
load("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffMeth/DiffMeth_Delta0.2.RData")
#load("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffExp/DiffExp_Comparative_analysis.RData")
rm(list= ls()[!(ls() %in% c("BMIQ.Meth", "results.dms", "results.dms.hyper", "results.dms.hypo","Tumor.BMIQ", "Normal.BMIQ"))])

######################### CpGs Annotation file from Hui Shen data ####################################
hm450 <- read.table("F:/OneDrive - University of Nebraska Medical Center/Human450K-Extra/HuiShen/April 2018/hm450.hg38.manifest.gencode.v22.tsv", sep="\t", header=TRUE, stringsAsFactors = FALSE, check.names = FALSE, row.names=5)
hm450 <- hm450[-c(1:5120),]

symbol <- strsplit(hm450$geneNames, split = ";")
type <- strsplit(hm450$transcriptTypes, split = ";")
transcript <- strsplit(hm450$transcriptIDs, split = ";")
positionTSS <- strsplit(hm450$distToTSS, split = ";")
feature <- strsplit(hm450$CGIposition, split = ";")
JHU.annotation <- data.frame(ProbeID = rep(rownames(hm450), sapply(positionTSS, length)), Chromosome=rep(hm450$CpG_chrm, sapply(positionTSS, length)),Strat=rep(hm450$CpG_beg, sapply(positionTSS, length)),End=rep(hm450$CpG_end, sapply(positionTSS, length)),Gene_Symbol = unlist(symbol),Transcript_ID = unlist(transcript), Gene_Type = unlist(type), Position_to_TSS=unlist(positionTSS), Feature_type = rep(hm450$CGIposition, sapply(positionTSS, length)))
JHU.annotation <- JHU.annotation[!is.na(JHU.annotation$Gene_Symbol),] ## Remove row which have <NA> in gene symbols
#tmp<- JHU.annotation[!is.na(JHU.annotation$Gene_Symbol),]
JHU.annotation$Position_to_TSS <-  varhandle::unfactor(JHU.annotation$Position_to_TSS)
JHU.annotation.TSS <-  JHU.annotation[((JHU.annotation$Position_to_TSS >= -1000) & (JHU.annotation$Position_to_TSS <= 500)),]
results.dms.hyper.TSS <- results.dms.hyper[rownames(results.dms.hyper)%in%JHU.annotation.TSS$ProbeID,]

rm(symbol, type, transcript, positionTSS, feature)

BMIQ.Meth.Hyper.TSS <- BMIQ.Meth[rownames(results.dms.hyper.TSS),]
# cancerID <- grep("-01A", colnames(BMIQ.Meth.Hyper.TSS))
# normalID <- grep("-11A", colnames(BMIQ.Meth.Hyper.TSS))
# BMIQ.Meth.Hyper.TSS <- BMIQ.Meth.Hyper.TSS[, cancerID]
clin <- GDCquery_clinic("TCGA-CHOL","clinical")

# Create directory (if it doesn't exist)
if (!dir.exists("Survival_Plot_DM_Hyper")) dir.create("Survival_Plot_DM_Hyper")   ### change directory name here

### R function for the survival analysis
plot_surv <- function (clinical_patient, dataMETH, Genelist, Survresult = FALSE, ThreshTop = 0.5, ThreshDown = 0.3, PercentUp = 25, PercentDown = 25, p.cut = 0.05) 
{
  
  if (!is.null(Genelist)) {
    Genelist <- intersect(rownames(dataMETH), Genelist)
    group1 <- colnames(dataMETH[,grep("-11A", colnames(dataMETH))])
    group2 <- colnames(dataMETH[,grep("-01A", colnames(dataMETH))])
    dataNormal <- dataMETH[Genelist, group1, drop = FALSE]
    dataCancer <- dataMETH[Genelist, group2, drop = FALSE]
    colnames(dataCancer) <- substr(colnames(dataCancer), 1, 12)
    numBeta.H <- rowSums(dataCancer > ThreshTop)
    numBeta.L <- rowSums(dataCancer < ThreshDown)
    #n1 <- rownames(Tumor.BMIQ[numBeta.H > round(ncol(Tumor.BMIQ)* (PercentUp/100)),])
    #n2 <- rownames(Tumor.BMIQ[numBeta.L > round(ncol(Tumor.BMIQ)* (PercentDown/100)),])
    dataCancer <- dataCancer[numBeta.H > round(ncol(dataCancer)* (PercentUp/100)) & numBeta.L > round(ncol(dataCancer)* (PercentDown/100)),]
    dataNormal <- dataNormal[rownames(dataCancer),]
  }
   if (is.null(Genelist)) {
    group1 <- colnames(dataMETH[,grep("-01A", colnames(dataMETH))])
    group2 <- colnames(dataMETH[,grep("-11A", colnames(dataMETH))])
    dataNormal <- dataMETH[, group1, drop = FALSE]
    dataCancer <- dataMETH[, group2, drop = FALSE]
    colnames(dataCancer) <- substr(colnames(dataCancer), 1, 12)
    numBeta.H <- rowSums(dataCancer > ThreshTop)
    numBeta.L <- rowSums(dataCancer < ThreshDown)
    #n1 <- rownames(Tumor.BMIQ[numBeta.H > round(ncol(Tumor.BMIQ)* (PercentUp/100)),])
    #n2 <- rownames(Tumor.BMIQ[numBeta.L > round(ncol(Tumor.BMIQ)* (PercentDown/100)),])
    dataCancer <- dataCancer[numBeta.H > round(ncol(dataCancer)* (PercentUp/100)) & numBeta.L > round(ncol(dataCancer)* (PercentDown/100)),]
    dataNormal <- dataNormal[rownames(dataCancer),]
  }
  
  cfu <- clinical_patient[clinical_patient[, "bcr_patient_barcode"] %in% substr(colnames(dataCancer), 1, 12), ]
  if ("days_to_last_followup" %in% colnames(cfu)) 
    colnames(cfu)[grep("days_to_last_followup", colnames(cfu))] <- "days_to_last_follow_up"
  cfu <- as.data.frame(subset(cfu, select = c("bcr_patient_barcode", "days_to_death", "days_to_last_follow_up", "vital_status")))
  if (length(grep("alive", cfu$vital_status, ignore.case = TRUE)) > 0) 
    cfu[grep("alive", cfu$vital_status, ignore.case = TRUE), "days_to_death"] <- "-Inf"
  if (length(grep("dead", cfu$vital_status, ignore.case = TRUE)) > 0) 
    cfu[grep("dead", cfu$vital_status, ignore.case = TRUE), "days_to_last_follow_up"] <- "-Inf"
  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]
  followUpLevel <- FALSE
  tabSurv_Matrix <- matrix(0, nrow(as.matrix(rownames(dataNormal))), 10)
  colnames(tabSurv_Matrix) <- c("mRNA", "pvalue", "Cancer Deaths", "Cancer Deaths with Top", "Cancer Deaths with Down", "Mean Tumor Top", "Mean Tumor Down", "Mean Normal", "Surv Mean Top Months", "Surv Mean Down Months")
  tabSurv_Matrix <- as.data.frame(tabSurv_Matrix)
  cfu$days_to_death <- as.numeric(as.character(cfu$days_to_death))
  cfu$days_to_last_follow_up <- as.numeric(as.character(cfu$days_to_last_follow_up))
  rownames(cfu) <- cfu[, "bcr_patient_barcode"]
  cfu <- cfu[!(is.na(cfu[, "days_to_last_follow_up"])), ]
  cfu <- cfu[!(is.na(cfu[, "days_to_death"])), ]
  cfu_complete <- cfu
  ngenes <- nrow(as.matrix(rownames(dataNormal)))
  for (i in 1:nrow(as.matrix(rownames(dataNormal)))) {
    cat(paste0((ngenes - i), "."))
    CpGselected <- as.matrix(rownames(dataNormal))[i]
    CpGselected_values <- dataCancer[rownames(dataCancer) == CpGselected, ]
    CpGselected_values_normal <- dataNormal[rownames(dataNormal) == CpGselected, ]
    if (all(CpGselected_values == 0)) 
      next
    tabSurv_Matrix[i, "mRNA"] <- CpGselected
    CpGselected_values_ordered <- sort(CpGselected_values, decreasing = TRUE)
    CpGselected_values_ordered_top <- as.numeric(ThreshTop)
    CpGselected_values_ordered_down <- as.numeric(ThreshDown)
    CpGselected_values_newvector <- CpGselected_values
    if (!is.na(CpGselected_values_ordered_top)) {
      numberOfSamples <- length(CpGselected_values_ordered)
      lastelementTOP <- max(which(CpGselected_values_ordered >= CpGselected_values_ordered_top))
      firstelementDOWN <- min(which(CpGselected_values_ordered <= CpGselected_values_ordered_down))
      samples_top_mRNA_selected <- names(CpGselected_values_ordered[1:lastelementTOP])
      samples_down_mRNA_selected <- names(CpGselected_values_ordered[firstelementDOWN:numberOfSamples])
      samples_UNCHANGED_mRNA_selected <- names(CpGselected_values_newvector[which((CpGselected_values_newvector) > 
                                                                                     CpGselected_values_ordered_down & CpGselected_values_newvector < 
                                                                                     CpGselected_values_ordered_top)])
      cfu_onlyTOP <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_top_mRNA_selected, ]
      cfu_onlyDOWN <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_down_mRNA_selected, ]
      cfu_onlyUNCHANGED <- cfu_complete[cfu_complete[, "bcr_patient_barcode"] %in% samples_UNCHANGED_mRNA_selected, ]
      cfu_ordered <- NULL
      cfu_ordered <- rbind(cfu_onlyTOP, cfu_onlyDOWN)
      cfu <- cfu_ordered
      ttime <- as.numeric(cfu[, "days_to_death"])
      sum(status <- ttime > 0)
      deads_complete <- sum(status <- ttime > 0)
      ttime_only_top <- cfu_onlyTOP[, "days_to_death"]
      deads_top <- sum(ttime_only_top > 0)
      if (dim(cfu_onlyDOWN)[1] >= 1) {
        ttime_only_down <- cfu_onlyDOWN[, "days_to_death"]
        deads_down <- sum(ttime_only_down > 0)
      }
      else {deads_down <- 0 }
      tabSurv_Matrix[i, "Cancer Deaths"] <- deads_complete
      tabSurv_Matrix[i, "Cancer Deaths with Top"] <- deads_top
      tabSurv_Matrix[i, "Cancer Deaths with Down"] <- deads_down
      tabSurv_Matrix[i, "Mean Normal"] <- mean(as.numeric(CpGselected_values_normal))
      dataCancer_onlyTop_sample <- dataCancer[, samples_top_mRNA_selected, drop = FALSE]
      dataCancer_onlyTop_sample_CpGselected <- dataCancer_onlyTop_sample[rownames(dataCancer_onlyTop_sample) == CpGselected, ]
      dataCancer_onlyDown_sample <- dataCancer[, samples_down_mRNA_selected, drop = FALSE]
      dataCancer_onlyDown_sample_CpGselected <- dataCancer_onlyDown_sample[rownames(dataCancer_onlyDown_sample) == CpGselected, ]
      tabSurv_Matrix[i, "Mean Tumor Top"] <- mean(as.numeric(dataCancer_onlyTop_sample_CpGselected))
      tabSurv_Matrix[i, "Mean Tumor Down"] <- mean(as.numeric(dataCancer_onlyDown_sample_CpGselected))
      ttime[!status] <- as.numeric(cfu[!status, "days_to_last_follow_up"])
      ttime[which(ttime == -Inf)] <- 0
      ttime <- ttime*0.032854884083862 ## conver days in months
      tabSurv_Matrix[i, "Surv Mean Top Months"] <- mean(ttime[1:nrow(cfu_onlyTOP)])
      tabSurv_Matrix[i, "Surv Mean Down Months"] <- mean(tail(ttime, nrow(cfu_onlyDOWN)))
      ttime1 <- ttime
      ttime <- Surv(ttime, status)
      rownames(ttime) <- rownames(cfu)
      legendHigh <- paste(CpGselected, "High")
      legendLow <- paste(CpGselected, "Low")
      legendHigh <- paste0(legendHigh, " (",nrow(cfu_onlyTOP), ")")
      legendLow <- paste0(legendLow, " (",nrow(cfu_onlyDOWN), ")")
      print(paste0("Now running analysis for the CpG :: ",CpGselected))
      tabSurv_pvalue <- tryCatch({
        tabSurv <- survdiff(ttime ~ c(rep("top", nrow(cfu_onlyTOP)), rep("down", nrow(cfu_onlyDOWN))))
        tabSurv_chis <- unlist(tabSurv)$chisq
        tabSurv_pvalue <- as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
      }, error = function(e) {
        return(Inf)
      })
      tabSurv_Matrix[i, "pvalue"] <- tabSurv_pvalue
      if (Survresult == TRUE) {
        
        titlePlot <- paste("Kaplan-Meier Survival analysis, pvalue=", round(tabSurv_pvalue, 3))
        plot(survfit(ttime ~ c(rep("high", nrow(cfu_onlyTOP)), rep("low", nrow(cfu_onlyDOWN)))), mark.time=TRUE, col = c("red4", "blue4"),  main = titlePlot, xlab = "Overall Survival Time", ylab = "Survival Probability", lwd = 1)
        legend("bottomleft", bty="n", lty=1, lwd=1.5, cex=0.8, legend = c(legendHigh, legendLow), col = c("red4","blue4"), text.col = c("red4","blue4"), pch = 15)
        #print(tabSurv)
        file_name = paste0("Survival_Plot_DM_Hyper/KM_Plot_", CpGselected, ".pdf")
        dev.print(pdf, file_name, width = 8, height = 8)
      }
    }
  }
  tabSurv_Matrix[tabSurv_Matrix == "-Inf"] <- 0
  tabSurvKM <- tabSurv_Matrix
  tabSurvKM <- tabSurvKM[tabSurvKM$mRNA != 0, ]
  tabSurvKM <- tabSurvKM[tabSurvKM$pvalue < p.cut, ]
  tabSurvKM <- tabSurvKM[!duplicated(tabSurvKM$mRNA), ]
  rownames(tabSurvKM) <- tabSurvKM$mRNA
  tabSurvKM <- tabSurvKM[, -1]
  tabSurvKM <- tabSurvKM[order(tabSurvKM$pvalue, decreasing = FALSE), ]
  return(tabSurvKM)
}
############################################################################

#Genelist <- rownames(BMIQ.Meth.Hyper.TSS)
results <- plot_surv(clin, dataMETH = BMIQ.Meth.Hyper.TSS, Genelist = rownames(BMIQ.Meth.Hyper.TSS), ThreshTop = 0.6, ThreshDown = 0.4, Survresult = TRUE, PercentUp = 25, PercentDown = 25, p.cut = 0.05) # I already defines group1= normal and group2=cancer
#results <- results[results$`Cancer Deaths with Top` >4 & results$`Cancer Deaths with Down` >4,]
#junk <- dir(path="Survival_Plot_DM/", pattern="KM_Plot_")
#file.remove(junk)
unlink("Survival_Plot_DM_Hyper/KM_Plot_*") ## Remove all Survival plot files
results.DM <- plot_surv(clin, dataMETH = BMIQ.Meth.Hyper.TSS, Genelist = rownames(results), ThreshTop = 0.6, ThreshDown = 0.4, Survresult = TRUE, PercentUp = 25, PercentDown = 25, p.cut = 0.05) # I already defines group1= normal and group2=cancer
#results <- results[results$`Cancer Deaths with Top` >4 & results$`Cancer Deaths with Down` >4,]

save.image("Survival_analysis_Promoter_Methylation_Hyper.RData")
