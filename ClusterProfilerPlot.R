library(clusterProfiler)
library(org.Hs.eg.db)
detach("package:dplyr", unload=TRUE)
gene <- DiffExp_edger_limma$geneID
geneList <- DiffExp_edger_limma$logFC
names(geneList) <- DiffExp_edger_limma$geneID
geneList <- sort(geneList, decreasing = TRUE)
ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "BP",
               level    = 3,
               readable = TRUE)

ego <- enrichGO(gene          = gene,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01,
                readable      = TRUE)
dotplot(ego)
barplot(ego, showCategory=10)

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
               minGSSize    = 5,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
summary(kk2)
barplot(kk, showCategory=10)

hsa00980 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa00980",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1), low = 4, mid = 8, high = 2, new.signature=FALSE,  res=600, key.pos="bottomleft")
##  high, low, mid can be use for color of plot