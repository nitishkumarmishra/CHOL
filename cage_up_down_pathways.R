library(dplyr)
library(pathview)
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
foldchange  <- DiffExp_edger_limma$logFC
names(foldchange) <- DiffExp_edger_limma$geneID
head(foldchange)
keggres = gage(foldchange, gsets=kegg.sets.hs, same.dir=TRUE)
lapply(keggres, head)
lapply(keggres, head, 10) # display top 10 pathways by up and downregulated genes

keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways
keggresids = substr(keggrespathways, start=1, stop=8)
detach("package:dplyr", unload=TRUE)## We need to detach dplyr, if loaded it's create proble in next step
plot_pathway = function(pid) pathview(gene.data=foldchange, pathway.id=pid, species="hsa", new.signature=FALSE)
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchange, pathway.id=pid, species="hsa"))

library(dplyr)
keggrespathways1 = data.frame(id=rownames(keggres$less), keggres$less) %>% 
  tbl_df() %>% 
  filter(row_number()<=5) %>% 
  .$id %>% 
  as.character()
keggrespathways1
keggresids1 = substr(keggrespathways1, start=1, stop=8)
detach("package:dplyr", unload=TRUE)
#plot_pathway = function(pid) pathview(gene.data=foldchange, pathway.id=pid, species="hsa", new.signature=FALSE)
tmp = sapply(keggresids1, function(pid) pathview(gene.data=foldchange, pathway.id=pid, species="hsa"))

# Manual blue and red color for plot
tmp = sapply(keggresids1, function(pid) pathview(gene.data=foldchange, pathway.id=pid, species="hsa", low = 4, mid = 8, high = 2, new.signature=FALSE,  res=600, key.pos="bottomleft"))

data(go.sets.hs)
data(go.subs.hs)
gobpsets = go.sets.hs[go.subs.hs$BP]

gobpres = gage(foldchange, gsets=gobpsets, same.dir=TRUE)

lapply(gobpres, head)

