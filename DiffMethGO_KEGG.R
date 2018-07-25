library(limma)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
DiffMethDeltaBeta.0.2 <- read.csv(file = "DiffMethDelta.0.2.txt", header = TRUE, sep = "\t", row.names = 1)
sigCPgS <- rownames(DiffMethDeltaBeta.0.2)
all <- rownames(BMIQ.Meth)
## For enrichment analysis I am using BMIQ.Meth probes as background rather than all CpGs
gstGO <- gometh(sig.cpg = sigCPgS, all.cpg = all, collection = "GO", plot.bias = TRUE)
gstKEGG <- gometh(sig.cpg = sigCPgS, all.cpg = all, collection = "KEGG", plot.bias = TRUE)

#gsaKEGG <- gsameth(sig.cpg = sigCPgS, all.cpg = all, collection="C:/Users/nitish.mishra/Desktop/MSigDB/c2.cp.kegg.v5.1.entrez.gmt", plot.bias = TRUE)
#gsaGO <- gsaGO <- gsameth(sig.cpg = sigCPgS, all.cpg = all, collection="C:/Users/nitish.mishra/Desktop/MSigDB/c1.all.v5.1.entrez.gmt", plot.bias = TRUE)
top20GO <- topGO(gstGO, number = 20, ontology = "BP")
top20KEGG <- topKEGG(gstKEGG, number = 20)


DiffMethDeltaBeta.0.3 <- read.csv(file = "DiffMethDelta.0.3.txt", header = TRUE, sep = "\t", row.names = 1)
sigCPgS_0.3 <- rownames(DiffMethDeltaBeta.0.3)
gstGO_0.3 <- gometh(sig.cpg = sigCPgS_0.3, all.cpg = all, collection = "GO", plot.bias = TRUE)
gstKEGG_0.3 <- gometh(sig.cpg = sigCPgS_0.3, all.cpg = all, collection = "KEGG", plot.bias = TRUE)
top20GO_0.3 <- topGO(gstGO_0.3, number = 20, ontology = "BP")
top20KEGG_0.3 <- topKEGG(gstKEGG_0.3, number = 20)


save.image(file = "DiffMethGO_KEGG.RData")
