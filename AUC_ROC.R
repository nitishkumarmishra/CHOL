require(pROC)
require(arsenal)

t.Meth <- t(Meth) ## Meth is methylation matrix
t.Meth <- as.data.frame(t.Meth)
t.Meth$class <- as.numeric(c(rep("1", 18),rep("2", 18))) ## Here I have 18 disease and and 18 normal samples in Meth.


datalist = list()
n <- ncol(t.Meth)-1
for(i in 1:n)
{
  fit <- glm(class ~ t.Meth[,i], data=t.Meth, family=binomial) ## This line will calculate AUC for all probe
  ## If you want to calculate specific probe e.g cg27384352 you have to use below command. Otherwise you can provide list of CpG's name also in this code.
  #fit <- glm(class ~ cg27384352, data=t.Meth, family=binomial)
  t.Meth$class <- as.integer(c(rep("0", 18),rep("1", 18)))
  tmp <- data.frame(summary(fit)$coef)
  tmp$OR <- round(exp(tmp[,1]),3)
  tmp$lower.CI <- round(exp(tmp[,1] +0.674* tmp[,2]),3)#75th percentile
  tmp$upper.CI <- round(exp(tmp[,1] +1.645* tmp[,2]),3)#95th percentile
  names(tmp)[4] <- 'P-value'
  pred <- predict(fit, type='response')
  #If plot=TRUE and smooth=TRUE it will draw smooth plot
  #tmp1 <- pROC::roc(t.Meth$class[!is.na(t.Meth$class)]~ pred, plot=TRUE, smooth=TRUE, percent=FALSE)### This line will draw ROC plot for each probe in present working directory
  tmp1 <- pROC::roc(t.Meth$class[!is.na(t.Meth$class)]~ pred, plot=FALSE, percent=FALSE) ## I am not using plot=FALSE, so it don't make ROC plot. 
  #May be you can select top CpG's and make list og CpG's and re-run this code nad make plot= TRUE. It will save time and also don't create to much ROC plot in the current folder.
  tmp1$auc[[1]]
  tmp$auc <- tmp1$auc[[1]]
  tt <- tmp[1,] ## This is the final results 
  rownames(tt) <- colnames(t.Meth[,i])
  datalist[[i]] <- tt
  
}
big_data = do.call(rbind, datalist)
rownames(big_data) <- colnames(t.Meth[,1:n])
big_data <- format(round(big_data, 3), nsmall = 2)
#rownames(big_data) <- colnames(t.Meth[,1:n]) ## if it don't display rownames, in loop may be
write.csv(big_data, file = "big_data_AUC.txt")