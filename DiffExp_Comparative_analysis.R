##aradegar is the list of genes with loogFC and adj.pvalue
##RPKM for a given GeneX is calculated by: (raw read counts * 10^9) / (total reads * length of GeneX).
#Input file for selecting top 1500 genes = *.uncv2.mRNAseq_RSEM_normalized_log2.txt (quantile normalized RSEM with log2 transformed) / *.mRNAseq_RPKM_log2.txt (RPKM value with log2 transformed) from mRNAseq_preprocess pipeline
library(pacman)
p_load('edgeR', 'DESeq2', 'limma')
geneExp <- read.csv("CHOL.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.data.txt", header = TRUE, sep = "\t", stringsAsFactors=FALSE, row.names = 1, check.names = FALSE)
geneExp <- geneExp[-c(1),which(geneExp[1,]=="raw_count")]
rownames(geneExp) <- gsub(pattern = "SLC35E2\\|728661", "SLC35E2B\\|728661", x = rownames(geneExp))
geneExp <- round(data.matrix(geneExp),0)
geneExp <- geneExp[-grep("\\?", rownames(geneExp)),] ## Remove genes which HGNC name is known i.e. ?
cancerID <- grep("01A", colnames(geneExp))
normalID <- grep("11A", colnames(geneExp))
cancerExp <- geneExp[,cancerID]
cnts <- cancerExp[rowSums(cancerExp==0)< ncol(cancerExp)*0.2,] ## Remove all gene which have 25% zero's
keep <- rowSums(cpm(cnts)>1) >= ncol(cnts)*0.25 #### Select only genes which have have CPM > 1 for >=50% samples
cnts <- cnts[keep,]
gene <- rownames(cnts)
cnts <- cbind(geneExp[gene,cancerID], geneExp[gene,normalID])
#Liver_genes <- read.csv(file = "Liver_Genes.txt", header = TRUE, sep = "\t")
#rownames(Liver_genes) <- Liver_genes$Genes

Liver <- read.csv(file = "Liver_Genes.txt", header = TRUE, sep = "\t")
Liver_ProteinAtlas <- read.csv(file = "ProteinAtlas_tissue_specificity_rna_liver_elevated.tab", header = TRUE, sep = "\t")
Liver_ProteinAtlas_Uniq <- Liver_ProteinAtlas[Liver_ProteinAtlas$RNA.tissue.category!='Group enriched',]
Liver_genes <- union(Liver$Genes, Liver_ProteinAtlas_Uniq$Gene)

#cnts.liver <- subset(cnts, !(rownames(cnts)%in%rownames(Liver_genes)))
#cnts.liver <- subset(cnts, !sapply(strsplit(rownames(cnts),"\\|"),'[[',1)%in%rownames(Liver_genes))
cnts.liver <- subset(cnts, !sapply(strsplit(rownames(cnts),"\\|"),'[[',1)%in%Liver_genes)
#####################Limma+Voom  #############
design <- model.matrix(~0 + factor(c(rep(2, length(cancerID)), rep(1, length(normalID)))))
colnames(design) <- c("Normal", "Tumor")
cont.matrix <- makeContrasts("Tumor-Normal", levels = design)
factors <- factor(c(rep("Tumor", length(cancerID)), rep("Normal", length(normalID))))
cnts.dgelist <- DGEList(cnts.liver, group=factors)
cnf <- calcNormFactors(cnts.dgelist, method = "TMM")
#cnf <- calcNormFactors(cnts.liver, method = "TMM") ## TMM normalization 
#v <- voom(cnts.liver, design, plot = TRUE, lib.size=colSums(cnts.liver) * cnf) ##Transform count data to log2-counts per million (logCPM)
v <- voom(cnf, design, plot = TRUE)
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
adj.method = "BH"; adj.pval = 0.01; raw.pval = 0.01; logFC = 1.5
aradeger <- topTable(fit2, adjust.method = adj.method, number = length(fit2), sort.by='logFC')
aradeger1 <- data.frame(aradeger[aradeger$adj.P.Val < adj.pval & aradeger$P.Value < raw.pval & abs(aradeger$logFC) >=logFC,])
diffExp.voom <- aradeger1
tmp <- decideTests(fit2, p.value = 0.01, lfc = 1.5)
summary(tmp)## We will get number of up and downregulated genes
diffExp.voom$geneID <- sapply(strsplit(rownames(diffExp.voom),"\\|"),'[[',2)
diffExp.voom$symbol<- sapply(strsplit(rownames(diffExp.voom),"\\|"),'[[',1)


###################### edgeR #####################
factors <- factor(c(rep("Tumor", length(cancerID)), rep("Normal", length(normalID))))
cnts.dgelist <- DGEList(cnts.liver, group=factors)
tmm <- calcNormFactors(cnts.dgelist, method = "TMM")
tmm <- estimateDisp(tmm) ## estimateDisp command do both CommonDisp and TagwiseDis
tmm.DE <- exactTest(tmm)
tmm.DE.top <- topTags(tmm.DE, n=nrow(tmm.DE), sort.by="logFC", p.value = 0.01)
topTags.DEG <- tmm.DE.top[tmm.DE.top$table$FDR < 0.01 & abs(tmm.DE.top$table$logFC) >=1.5,]
diffExp.edgeR <- topTags.DEG$table
diffExp.edgeR$geneID <- sapply(strsplit(rownames(diffExp.edgeR),"\\|"),'[[',2)
diffExp.edgeR$symbol<- sapply(strsplit(rownames(diffExp.edgeR),"\\|"),'[[',1)
###################### Deseq2 #####################
cond <- factor(c(rep(2, length(cancerID)), rep(1, length(normalID))))
dds <- DESeqDataSetFromMatrix(cnts.liver, DataFrame(cond), ~ cond)
dds <- DESeq(dds)
p_unload(GenomicDataCommons)#unload GenomicDataCommons before results, because DESeq2 and GenomicDataCommons both have results command
#resmiRNA <- results(dds, lfcThreshold = 2, pAdjustMethod = "BH")
resmiRNA <- results(dds, pAdjustMethod = "BH")
resmiRNAOrdered <- resmiRNA[order(resmiRNA$log2FoldChange, decreasing = TRUE),]
resmiRNASig <- subset(resmiRNAOrdered, (padj < 0.01& abs(resmiRNAOrdered$log2FoldChange) >=1.5))
diffExp.DESeq2 <- as.data.frame(resmiRNASig)
diffExp.DESeq2$foldchange <- 2^diffExp.DESeq2$log2FoldChange
diffExp.DESeq2$geneID <- sapply(strsplit(rownames(diffExp.DESeq2),"\\|"),'[[',2)
diffExp.DESeq2$symbol<- sapply(strsplit(rownames(diffExp.DESeq2),"\\|"),'[[',1)
####################################################
list_of_data = list(diffExp.edgeR , diffExp.DESeq2, diffExp.voom)
common_names = Reduce(intersect, lapply(list_of_data, row.names))
### common_names is list of gene names which is common in all three matrix
list_of_data = lapply(list_of_data, function(x) { x[row.names(x) %in% common_names,] })
lapply(list_of_data, nrow) ## to check the length of each data frame
#####################################################
save.image("DiffExp_Comparative_analysis.RData")
#####################################################
#Make input file for GSEA analysis
#Input <- -log10(aradeger$adj.P.Val)*sign(aradeger$logFC)
#Input <- cbind(sapply(strsplit(rownames(aradeger),"\\|"),'[[',2), Input)
Input <- -log10(resmiRNA$pvalue)*sign(resmiRNA$log2FoldChange)
Input <- cbind(sapply(strsplit(rownames(resmiRNA),"\\|"),'[[',1), Input)
### Load MsigDB data###############
load(file = "C://Users/nitish.mishra/Desktop/MSigDB/MSigDB5.1 RData/human_c2_v5p2.rdata")

y <- log2(cnts.liver+1)
rownames(y) <- sapply(strsplit(rownames(y),"\\|"),'[[',2)
idx <- ids2indices(Hs.c2, id=rownames(y))
cam <- camera(y, idx, design = design, contrast = cont.matrix, inter.gene.cor = 0.01)
cam <- cam[cam$PValue < 0.01 & cam$FDR < 0.05,]
cam <- cam[grep("KEGG|REACTOME|PID|BIOCARTA", rownames(cam)),]
options(digits = 3)
head(cam, 20)
####alternate option for KEGG, BIOCARTA, REACTOME, PID only ####
Hs.c2.new <- Hs.c2[grep("KEGG|REACTOME|PID|BIOCARTA", names(Hs.c2))]#### make file for only KEGG
idx <- ids2indices(Hs.c2.new, id=rownames(y))
cam <- camera(y, idx, design = design, contrast = cont.matrix, inter.gene.cor = 0.01)
cam <- cam[cam$PValue < 0.01 & cam$FDR < 0.05,]
head(cam)

################## barcode plot ######################
## Convert 
res <- tmm; rownames(res$counts) <- sapply(strsplit(rownames(cnts.liver),"\\|"),'[[',2)
res.DE <- exactTest(res)
barcodeplot(res.DE$table$logFC, index = idx[['KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION']], labels = c("Tumor", "Normal"), main="KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION")
barcodeplot(res.DE$table$logFC, index = idx[['CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_UP']], index2 = idx[['CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION_DN']],labels = c("Tumor", "Normal"), main="CHIANG_LIVER_CANCER_SUBCLASS_PROLIFERATION", alpha = 1)

## Function to remove genes with more than 50% zero's
# rem <- function(x){
#   x <- as.matrix(x)
#   x <- t(apply(x,1,as.numeric))
#   r <- as.numeric(apply(x,1,function(i) sum(i == 0)))
#   remove <- which(r >= dim(x)[2]*0.2)
#   return(remove)
# }
##########ClusterProfiler ##############################
library(clusterProfiler)
library(org.Hs.eg.db)
detach("package:dplyr", unload=TRUE)

diffExp.ClusterProfiler <- as.data.frame(resmiRNA)
diffExp.ClusterProfiler$geneID <- sapply(strsplit(rownames(diffExp.ClusterProfiler),"\\|"),'[[',2)
diffExp.ClusterProfiler$symbol<- sapply(strsplit(rownames(diffExp.ClusterProfiler),"\\|"),'[[',1)
diffExp.ClusterProfiler$weight <- -log10(diffExp.ClusterProfiler$padj)*sign(diffExp.ClusterProfiler$log2FoldChange)
#diffExp.ClusterProfiler$weight <- diffExp.ClusterProfiler$log2FoldChange
#diffExp.ClusterProfiler$weight <- scales::rescale(diffExp.ClusterProfiler$weight, to = c(-5, 5))
diffExp.ClusterProfiler <- diffExp.ClusterProfiler[order(diffExp.ClusterProfiler$weight, decreasing = TRUE),]

gene <- diffExp.DESeq2$geneID
geneList <- diffExp.ClusterProfiler$weight
names(geneList) <- diffExp.ClusterProfiler$geneID
#geneList <- sort(geneList, decreasing = TRUE)## Already sorted on basis of weight
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)
barplot(ego, drop=TRUE, showCategory=12)
bp <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)## Remove similar GO
barplot(bp, drop=TRUE, showCategory=12)

ego1 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              nPerm        = 1000,
              minGSSize    = 10,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
dotplot(ego)
barplot(ego, showCategory=10)


kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.01)
barplot(kk, showCategory=10)
dotplot(kk)

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               verbose      = FALSE)
summary(kk2)
head(kk2)
gseaplot(kk2, geneSetID = "hsa04110", title = "GSEA plot for Cell cycle (hsa04110)")
gseaplot(kk2, geneSetID = "hsa00040", title = "GSEA plot for Pentose and glucuronate interconversions (hsa04110)")

library(pathview)
geneList1 <- diffExp.ClusterProfiler$log2FoldChange
names(geneList1) <- diffExp.ClusterProfiler$geneID
pathview(gene.data  = geneList1,
         pathway.id = "hsa00020",
         species    = "hsa", 
         low = list(gene = "green", cpd = "blue"), 
         mid =list(gene = "gray", cpd = "gray"), 
         high = list(gene = "red", cpd = "yellow"),
         limit      = list(gene=max(abs(geneList1)), cpd=1))
file.remove("hsa00020.PNG", "hsa00020.XML")

avid <- enrichDAVID(gene = gene, idType = "ENTREZ_GENE_ID", listType = "Gene",
                     annotation = "KEGG_PATHWAY", david.user = "nitish.mishra@unmc.edu")## May be we need run two time, because email ID did't work first time

edo <- enrichDO(gene = gene, ont = "DO", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH")
edo1 <- gseDO(geneList = geneList, nPerm = 1000, minGSSize = 5, maxGSSize = 500, pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)
gseaplot(edo1, geneSetID = "DOID:1324")

### GSEA by using clusterProfiler #####
gmtfile <- system.file("extdata", "c5.mf.v5.2.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)

egmt <- enricher(gene, TERM2GENE=c5, universe =  names(geneList),
                 minGSSize = 10, 
                 maxGSSize = 500,
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(egmt)

egmt2 <- GSEA(geneList, TERM2GENE=c5, verbose=FALSE, nPerm = 1000,
              minGSSize = 10, 
              maxGSSize = 500,
              pAdjustMethod = "BH",
              pvalueCutoff  = 0.05)
head(egmt2)

gmtfile <- system.file("extdata", "c2.cp.kegg.v5.2.entrez.gmt", package="clusterProfiler")
c2 <- read.gmt(gmtfile)
kk3 <- GSEA(geneList, TERM2GENE=c2, verbose=FALSE, nPerm = 1000,
     minGSSize = 10, 
     maxGSSize = 500,
     pAdjustMethod = "BH",
     pvalueCutoff  = 0.05)
head(kk3)

save.image("CAMERA_DiffExp_analysis.RData")
############# Heatmap plot of top 50 DEGs ##############
logCPM <- cpm(cnts.dgelist, prior.count = 2, log = TRUE)
colnames(logCPM) <- substr(colnames(logCPM), 9, 15)
#o <- order(tmm.DE$table$PValue)
o <- order(resmiRNA$padj)
logCPM <- logCPM[o[1:50],]
logCPM <- t(scale(t(logCPM)))
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")
heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both")
# heatmap.2(logCPM, col=col.pan, Rowv=TRUE, scale="none", trace="none", dendrogram="both", 
#           cexRow=1, cexCol=1.4, density.info="none", margin=c(10,9), lhei=c(2,10), lwid=c(2,6))