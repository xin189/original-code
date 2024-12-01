library(dbplyr)
library(clusterProfiler)
library(AnnotationDbi)
library(AnnotationHub)
library(rtracklayer)

library(DOSE)
#library("org.Mm.eg.db")
args = commandArgs(T)
goat=loadDb("/scRNA/USER/zhouli/bin/Database/goat.OrgDb")
data=read.table(args[1],h=F,row.names=1,skip =1)
up=data[data$V3 > 0,]
down=data[data$V3 < 0,]
sp=list(up=up,down=down)
for(i in 1:length(names(sp))){
id<-bitr(row.names(sp[[i]]), fromType="SYMBOL", toType="ENTREZID",OrgDb=goat)
ego<-enrichGO(OrgDb=goat, gene = id$ENTREZID,pvalueCutoff = 0.05,readable=TRUE,ont = "BP")
ekk<-enrichKEGG(gene=id$ENTREZID,keyType='kegg',organism='chx',pvalueCutoff=0.05) #http://www.genome.jp/kegg/catalog/org_list.html
plotdata<-ego@result
plotdata$GeneRatio<-sapply(plotdata$GeneRatio,function(x) eval(parse(text = x)))
write.table(plotdata,paste0(args[2],".",names(sp)[i],".ego.txt"),sep="\t",row.names=FALSE,quote=FALSE)
plotdata<-ekk@result
plotdata$geneID <- sapply(strsplit(as.character(plotdata$geneID), "/"), function(x) {
          entrez_ids <- as.integer(x)
          entrez_ids <- as.character(entrez_ids)
          symbols <- select(goat, keys = entrez_ids, columns = "SYMBOL", keytype = "ENTREZID")$SYMBOL
          return(paste(symbols, collapse = "/"))
        })
plotdata$GeneRatio<-sapply(plotdata$GeneRatio,function(x) eval(parse(text = x)))
write.table(plotdata,paste0(args[2],".",names(sp)[i],".ekk.txt"),sep="\t",row.names=FALSE,quote=FALSE)
pdf(paste0(args[2],".",names(sp)[i],".GO_KEGG.enrich.pdf"),h=8,w=8)
print(barplot(ego,showCategory=11,drop=T))
print(dotplot(ego,showCategory=10))
print(barplot(ekk,showCategory=10,drop=T))
print(dotplot(ekk,showCategory=10))
eGoBP <- enrichGO(gene = id$ENTREZID,
        OrgDb = goat, 
        pvalueCutoff =0.01, 
        qvalueCutoff = 0.01,
        ont="BP", #BP\MF\CC
        readable =T)
print(plotGOgraph(eGoBP))
dev.off()
}
