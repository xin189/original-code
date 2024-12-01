library(monocle)
library(dplyr)
library(Seurat)
args = commandArgs(T)

combined_sub<- readRDS(file=args[1])
#combined_sub = combined_sub[,sample(colnames(combined_sub),20000)]
#combined_sub = subset(combined_sub,seurat_clusters %in% c('1','2'))
setwd(args[2])
pdf("subset_celltype.pdf")
DimPlot(object = combined_sub,reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 3)+theme(legend.position = 'bottom')+theme(legend.text = element_text(size = 11))
dev.off()
Idents(combined_sub)=combined_sub@meta.data$seurat_clusters
combined_sub@meta.data$seurat_clusters<-Idents(combined_sub)  ######
data<- as(as.matrix(combined_sub@assays$RNA@counts),'sparseMatrix')  #for integrated data
write.table(combined_sub@assays$RNA@counts,file="sub_count.csv",quote=F)
write.table(combined_sub@meta.data["stim","seurat_clusters"],file="sub_cluster.csv",quote=F)


pd<-combined_sub@meta.data
fd<-data.frame(gene_short_name = row.names(data),row.names=row.names(data))
pd <- new("AnnotatedDataFrame", data = pd)
fd <- new("AnnotatedDataFrame", data = fd)
my_cds <- newCellDataSet(data,phenoData = pd, featureData = fd)

#my_cds <- newCellDataSet(as.matrix(gbm),phenoData = pd,featureData = fd,lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())
pData(my_cds)$barcode<-colnames(combined_sub)

my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)
my_cds <- detectGenes(my_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= 3))
my_cds_subset <- my_cds[expressed_genes,]
head(pData(my_cds_subset))

disp_table <- dispersionTable(my_cds_subset) 
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.01)
my_cds_subset <- setOrderingFilter(my_cds_subset, unsup_clustering_genes$gene_id)

# standardise to Z-distribution(标准正态分布)
x <- pData(my_cds_subset)$num_genes_expressed
x_1 <- (x - mean(x)) / sd(x)
library(ggplot2)
library(cowplot)
df <- data.frame(x = x_1)
pdf(file="expressed_genes_distrebute.pdf") 
ggplot(df, aes(x)) +  geom_histogram(bins = 50) +  geom_vline(xintercept = c(-2, 2), linetype = "dotted", color = 'red')
dev.off()

pData(my_cds_subset)$UMI <- Matrix::colSums(exprs(my_cds_subset))
pdf(file="expressed_genes_UMI.pdf") 
ggplot(pData(my_cds_subset), aes(num_genes_expressed, UMI)) + geom_point()
dev.off()

#Cluster
#expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= 3))
#my_cds_subset <- my_cds[expressed_genes,]
#my_cds_subset <- detectGenes(my_cds_subset, min_expr = 0.1)
#head(pData(my_cds_subset))

#Constructing Single Cell Trajectories
clustering_DEG_genes <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = '~seurat_clusters',cores = 2)
#my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000] 
my_ordering_genes<-row.names(subset(clustering_DEG_genes, qval < 0.01))
my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = my_ordering_genes)
pdf(file="Monocle_ordering_genes.pdf")
plot_ordering_genes(my_cds_subset)
dev.off()
my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree',auto_param_selection = F)
my_cds_subset <- orderCells(my_cds_subset)

#GM_state <- function(cds){
#  if (length(unique(pData(cds)$State)) > 1){
#    T0_counts <- table(pData(cds)$State, pData(cds)$seurat_clusters)[,"OPC"]
#    return(as.numeric(names(T0_counts)[which
#                                       (T0_counts == max(T0_counts))]))
#} else {
#  return (1)
#}}
#
#my_cds_subset <- orderCells(my_cds_subset,root_state = GM_state(my_cds_subset),num_paths =5)

write.table(pData(my_cds_subset)[,c("barcode","num_genes_expressed","UMI","seurat_clusters","State","Pseudotime")],row.names = F,file="Monocle_Cell_Summary_baseonCluster.xls",sep="\t",quote = F)




pdf(file="Monocle_Cell_Trajectories_seurat_clusters.pdf")
plot_cell_trajectory(my_cds_subset, color_by = "seurat_clusters",show_branch_points = F)
#plot_cell_trajectory(my_cds_subset, color_by = "as.factor(seurat_clusters)")
dev.off()

pdf("Cell_Trajectories_seurat_clusters_facet.pdf")
plot_cell_trajectory(my_cds_subset, color_by = "seurat_clusters")+facet_wrap(~seurat_clusters, nrow = 1)
dev.off()

pdf("Pseudotime.pdf")
plot_cell_trajectory(my_cds_subset, color_by = "Pseudotime",show_branch_points = F)
dev.off()

pdf(file="Monocle_Cell_Trajectories_State.pdf") 
plot_cell_trajectory(my_cds_subset, color_by = "State",show_branch_points =F)
dev.off()

saveRDS(my_cds_subset,file = "monocle2_subset.rds")
q()
## diff genes from diff clusters
diff_cluster<-differentialGeneTest(my_cds_subset,fullModelFormulaStr = '~seurat_clusters',cores = 8)
order_genes<-row.names(subset(clustering_DEG_genes, qval < 0.01))
write.table(clustering_DEG_genes[order_genes,], file='cluster_de_filter.xls', sep='\t',quote = F)
pdf("cluster_pseudotime_heatmap.pdf")
plot_pseudotime_heatmap(my_cds_subset[order_genes,],num_clusters = 3,cores = 1,show_rownames = T)
dev.off()
pdf("cluster_in_pseudotime.pdf")
clustering_DEG_genes %>% subset(qval < 0.01) %>% arrange(qval) %>% head() -> top6_gene
top6_gene<-rownames(top6_gene)
plot_genes_in_pseudotime(my_cds_subset[top6_gene,], color_by = "seurat_clusters")
dev.off()

##### diff genes from diff State
diff_State<- differentialGeneTest(my_cds_subset, fullModelFormulaStr = "~State", cores = 8)
#sig_gene_names <- row.names(subset(diff_State, qval < 0.01))
diff_State %>% subset(qval < 0.01) %>% arrange(qval) %>% head(20) -> sig_gene_names
write.table(sig_gene_names, file='my_state_de_head20.xls', sep='\t',quote = F)
sig_gene_names<- rownames(sig_gene_names)

#sig_gene_names<-head(sig_gene_names)
pdf("State_pseudotime_heatmap.pdf")
pseudotime_heatmap<-plot_pseudotime_heatmap(my_cds_subset[sig_gene_names,],num_clusters = 3,cores = 1,show_rownames = T)
head(pseudotime_heatmap)
dev.off()
pdf("State_in_pseudotime.pdf")
plot_genes_in_pseudotime(my_cds_subset[sig_gene_names,], color_by = "State")
dev.off()


#Finding Genes that Change as a Function of Pseudotime
my_pseudotime_de <- differentialGeneTest(my_cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 8)
my_pseudotime_de %>% arrange(qval) %>% head()-> head_pseudotime_gene
head_pseudotime_gene <- rownames(head_pseudotime_gene)
pdf(file="Monocle_Top6_significant_Diff_genes_in_pseudotime.pdf")
plot_genes_in_pseudotime(my_cds_subset[head_pseudotime_gene,])
dev.off()
write.table((my_pseudotime_de %>% arrange(qval))[,c("gene_short_name","pval","qval")],row.names = F,col.names = c("GeneName","pval","qval"),file="Monocle_Cluster_Pseudotime.GeneDiffExp.xls",sep="\t",quote = F)
write.table((my_pseudotime_de %>% arrange(qval) %>% head(50))[,c("gene_short_name","pval","qval")],row.names = F,col.names = c("GeneName","pval","qval"),file="Monocle_Cluster_Pseudotime.GeneDiffExpFilter.xls",sep="\t",quote = F)

my_pseudotime_de %>% arrange(qval) %>% head(50) -> gene_to_cluster
gene_to_cluster <- rownames(gene_to_cluster)
pdf(file="Top50_significant_genes_in_pseudotime_heatmap.pdf") 
my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,], num_clusters = 3, cores = 8, show_rownames = TRUE, return_heatmap = TRUE)
dev.off()

#Analyzing Branches in Single-Cell Trajectories
#for (i in unique(Branchesdata$branch_point_idx)){
for (i in 1:length(my_cds_subset@auxOrderingData$DDRTree$branch_points)){
	BEAM_res <- BEAM(my_cds_subset, branch_point = i, cores = 8)
	BEAM_res <- BEAM_res[order(BEAM_res$qval),]
	res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
	BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
	pdfname <- paste("Monocle_Top100_significant_genes_in_branched",i,"_heatmap.pdf",sep='')
	pdf(file=pdfname)
	my_branched_heatmap <- plot_genes_branched_heatmap(my_cds_subset[row.names(head(BEAM_res, 100)),],branch_point = i,num_clusters = 4,cores = 8,use_gene_short_name = TRUE,show_rownames = TRUE,return_heatmap = TRUE)
	dev.off()
	filename1<-paste("Monocle_Cluster_Branches",i,".GeneDiffExp.xls",sep='')
	filename2<-paste("Monocle_Cluster_Branches",i,".GeneDiffExpFilter.xls",sep='')
	write.table(res,file=filename1,row.names = F,col.names = c("GeneName","pval","qval"),sep="\t",quote = F)
	write.table(head(res,100),file=filename2,row.names = F,col.names = c("GeneName","pval","qval"),sep="\t",quote = F)
}
#monocle对象提取子集
#monocle_sub<-subset(my_cds_subset,pData(my_cds_subset)$stim=='test')
