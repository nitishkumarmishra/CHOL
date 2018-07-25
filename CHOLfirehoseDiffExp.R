library(DESeq2)
##### R code for analysis of TCGA CHOL data from Firehose ###########
geneExp <- read.csv("CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt", header = TRUE, sep = "\t", row.names = 1)
#geneExp <- geneExp[-1,]
colnames(geneExp) <- gsub("\\.", "-", colnames(geneExp))
miRNAExp <- read.csv("CHOL.mirnaseq__illuminahiseq_mirnaseq__bcgsc_ca__Level_3__miR_gene_expression__data.data.txt", header = TRUE, sep = "\t", row.names = 1)
miRNAExp <- miRNAExp[,seq(1,135,3)]### read colum after 3, from 1 to 135
#miRNAExp <- miRNAExp[-1,]
colnames(miRNAExp) <- gsub("\\.", "-", colnames(miRNAExp))
####################################################################
colnames(geneExp) <- substr(colnames(geneExp), 1, 16)
colnames(miRNAExp) <- substr(colnames(miRNAExp), 1, 16)
save(geneExp, file = "CHOL-firehoseExp.Rda")
save(miRNAExp, file = "CHOL-firehoseMiRNA.Rda")
cancerID <- grep("01A", colnames(miRNAExp))
normalID <- grep("11A", colnames(miRNAExp))
####################################################################
cnts <- cbind(miRNAExp[,cancerID], miRNAExp[,normalID])
cnts <- data.matrix(cnts)
#cnts <- cnts[rowSums(cnts) > 0, ] ### remove row with all zero
### Remove sample which have more than 20% zero, similarly sample which have more than 20% zero
cnts <- cnts[apply(cnts,1,function(x) sum(x==0))<ncol(cnts)*0.8,]
cnts <- cnts[apply(cnts,2,function(x) sum(x==0))<nrow(cnts)*0.8,]
cond <- factor(c(rep(2, length(cancerID)), rep(1, length(normalID))))
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
dds <- DESeq(dds)
resmiRNA <- results(dds, lfcThreshold = 1, pAdjustMethod = "BH")
resmiRNAOrdered <- resmiRNA[order(resmiRNA$padj),]
resmiRNASig <- subset(resmiRNAOrdered, padj <= 0.05)
tmp <- as.data.frame(resmiRNASig)
tmp$foldchange <- 2^tmp$log2FoldChange
write.csv(as.data.frame(tmp),file="miRNA_results.csv")
##write.csv(as.data.frame(resSig), file="miRNA_results.csv")
# counts.norm = data.frame(counts(dds, normalized=TRUE))
# cond1 = counts.norm[groups$cond == 1]
# cond2 = counts.norm[groups$cond == 2]
# cond1.means = apply(cond1,1,mean)
# cond2.means = apply(cond2,1,mean)

# # munge data and output
# dat = data.frame(rownames(res),
#                  round(res$log2FoldChange,digits=3),
#                  signif(res$padj,digits=3),
#                  round(cond1.means,digits=3), 
#                  round(cond2.means,digits=3), 
#                  round((cond1.means/cond2.means), digits = 3))
# names(dat) <- c('GeneID','log2foldchange','significance','Tumor','Normal', 'FoldChange')
# write.csv(dat, file="outFile.csv",sep="\t",row.names=F)
#######################################################################
cancerID <- grep("01A", colnames(geneExp))
normalID <- grep("11A", colnames(geneExp))
cnts <- cbind(geneExp[,cancerID], geneExp[,normalID])
cnts <- data.matrix(cnts)
#cnts <- cnts[rowSums(cnts) > 0, ]
### Remove sample which have more than 20% zero, similarly sample which have more than 20% zero
cnts <- cnts[apply(cnts,1,function(x) sum(x==0))<ncol(cnts)*0.8,]
cnts <- cnts[apply(cnts,2,function(x) sum(x==0))<nrow(cnts)*0.8,]
cnts <- round(cnts)
cond <- factor(c(rep(2, length(cancerID)), rep(1, length(normalID))))
dds <- DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
dds <- DESeq(dds)
res <- results(dds, lfcThreshold = 1, pAdjustMethod = "BH")
resOrdered <- res[order(res$padj),]
resSig <- subset(resOrdered, padj <= 0.01)
tmp <- as.data.frame(resSig)
tmp$geneID <- sapply(strsplit(rownames(tmp),"\\|"),'[[',2)
rownames(tmp) <- sapply(strsplit(rownames(tmp),"\\|"),'[[',1)
tmp$foldchange <- 2^tmp$log2FoldChange
write.csv(tmp,file="mRNA_results.0.01.csv")

