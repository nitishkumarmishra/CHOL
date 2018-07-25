setwd("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffExp")
load("F:/OneDrive - University of Nebraska Medical Center/eMap/RTCGAToolbox/TCGABiolink/dataTGCT/CHOL/DiffExp/DiffExp_Comparative_analysis.RData")
rm(list= ls()[!(ls() %in% c("rna", "common_names", "geneExp"))])
clin <- GDCquery_clinic("TCGA-CHOL","clinical")
rna <- log2(cnts.liver+1)
cancerID <- grep("-01A", colnames(rna))
normalID <- grep("-11A", colnames(rna))
rna <- cbind(rna[,cancerID], rna[,normalID])

rownames(rna) <- sapply(strsplit(rownames(rna),"\\|"),'[[',1)
common_id <- sapply(strsplit(common_names,"\\|"),'[[',2)
common_names <- sapply(strsplit(common_names,"\\|"),'[[',1)
common_list <- as.data.frame(cbind(common_names, common_id))
rownames(common_list) <- common_list$common_names

#dataCancer <- rna[common_names, cancerID, drop = FALSE]
#colnames(dataCancer) <- substr(colnames(dataCancer), 1, 12)

if (!dir.exists("Survival_Plot_DEG_1")) dir.create("Survival_Plot_DEG_1")   ### change directory name here
#load("../Survival_AUC_Final.RData")
plot_surv <- function (clinical_patient, dataGE, Genelist, Median = TRUE, Mean =TRUE) 
{
  Genelist <- intersect(rownames(dataGE), Genelist)
  group1 <- colnames(dataGE[,grep("-11A", colnames(dataGE))])
  group2 <- colnames(dataGE[,grep("-01A", colnames(dataGE))])
  dataNormal <- dataGE[Genelist, group1, drop = FALSE]
  dataCancer <- dataGE[Genelist, group2, drop = FALSE]
  colnames(dataCancer) <- substr(colnames(dataCancer), 1, 12)
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
  
  rownames(cfu) <- cfu$bcr_patient_barcode
  cfu$Days <- as.integer(ifelse(cfu$vital_status == "alive", cfu$days_to_last_follow_up, cfu$days_to_death))
  cfu$Time <- cfu$Days*0.032854884083862
  cfu$Status <- ifelse(cfu$vital_status=="alive", 0, 1)
  cfu_temp <- cfu
  #ngenes <- nrow(as.matrix(rownames(dataCancer)))
for(i in 1:nrow(common_Exp.AUC)){
  ngenes <- nrow(as.matrix(rownames(dataCancer)))
  cat(paste0((ngenes - i), "."))
  dataCancer <- dataCancer[rownames(dataCancer)%in% rownames(common_Exp.AUC),]
  mRNAselected <- as.matrix(rownames(dataCancer))[i]
  mRNAselected_values <- dataCancer[rownames(dataCancer) == mRNAselected, ]
  if(Mean==TRUE)
  {
     cfu_temp$Exp <- ifelse(mRNAselected_values > mean(mRNAselected_values), "High", "Low")
  }
 if(Median==TRUE)
 {
   cfu_temp$Exp <- ifelse(mRNAselected_values > median(mRNAselected_values), "High", "Low")
 }
  print(paste0("Now running analysis for the Gene :: ",mRNAselected, "::", tabSurv_pvalue))
  tabSurv <- survdiff(Surv(cfu_temp$Time, cfu_temp$Status) ~ ifelse(mRNAselected_values > mean(mRNAselected_values), "High", "Low"))
  tabSurv_chis <- unlist(tabSurv)$chisq
  tabSurv_pvalue <- as.numeric(1 - pchisq(abs(tabSurv$chisq), df = 1))
  print(tabSurv_pvalue)
  file_name = paste0("Survival_Plot_DEG_1/KM_Plot_", mRNAselected, ".pdf")
  if(tabSurv_pvalue <= 0.05)
  {
    
    #file_name = paste0("Survival_Plot_DEG_1/KM_Plot_", mRNAselected, ".pdf")
    #dev.print(pdf, file_name, width = 8, height = 8)
    #sfit <- survfit(Surv(cfu_temp$Time, cfu_temp$Status) ~ ifelse(mRNAselected_values > mean(mRNAselected_values), "High", "Low"), data=cfu_temp)
    sfit <- do.call(survfit, list(formula = Surv(Time, Status == 1) ~ Exp, data = cfu_temp))
    p <- ggsurvplot(sfit, pval = TRUE, risk.table = "abs_pct",  xlim = c(0, max(cfu$Time)+2), risk.table.height=0.15, risk.table.fontsize=4, pval.coord = c(0,0),legend.title= "Exp", pval.size = 5,  main="Kaplan-Meier Overall Survival Curves", ylab = "Probability of survival", palette=c("red4", "blue4"))
    #dev.print(pdf, file_name, width = 8, height = 8)
    pdf(file_name, width = 8, height = 8, onefile = FALSE)
    print(p)
    dev.off()
    #sfit <- do.call(survfit, list(formula = Surv(Time, Status) ~ Exp, data = cfu_temp))
    #ggsurvplot(sfit, pval = TRUE, risk.table = "abs_pct", risk.table.height=0.15, risk.table.fontsize=3, legend.title= "Exp", surv.median.line = "hv", pval.size = 5,  main="Kaplan-Meier Overall Survival Curves", ylab = "Probability of survival", palette=c("red4", "blue4"))
    #dev.print(pdf, file_name, width = 8, height = 8)
  }
  #if(tabSurv_pvalue > 0.05){next}
  }
}
 

plot_surv(clin, rna, Genelist = rownames(DiffExp_Surv_0.05), Median = TRUE)