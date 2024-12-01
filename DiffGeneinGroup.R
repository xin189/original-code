## visualize these changes in gene expression is with the split.by option 
library(tidyverse)
library(Seurat)
library(dplyr)
library(cowplot)
args <- commandArgs(T)
#top 6 cell

#load('/scRNA/USER/zhouli/project/tmp/ER047.wengshijun/result/process/vis_changes_gene/Merge1/old--vs--D0/combined.markers.Rdata')
combined=readRDS(args[1])
Idents(combined)='celltype'
cell_Type_name <- NA

if (1){
#if (is.na(cell_Type_name)){
	table_ident <- table(combined@active.ident)
	cell_Type_name = row.names(table_ident[table_ident>=3])

}else{
	table_ident <- table(combined@active.ident)
	cell_Type_name_all = row.names(table_ident[table_ident>=3])
	#<3 can not find marker
	cell_Type_name <- strsplit(cell_Type_name ,',')[[1]]
	cell_Type_name <- cell_Type_name[cell_Type_name %in% cell_Type_name_all]
	

}

#####
if(is.null(combined@meta.data$group)){
                combined_subset <- subset(combined , subset =stim %in% c(args[3],args[4]))
#        }else{
#                combined_subset <- subset(combined , subset =group %in% c("D0","old"))
#			combined_subset$stim<-combined_subset$group
        }
combined<-combined_subset

for ( i  in cell_Type_name ){
#combined.tmp <- combined

abc<-function(combined){
combined$celltype.stim <- paste(Idents(object = combined), combined$stim, sep = "_")
combined$celltype <- Idents(object = combined)
Idents(object = combined) <- "celltype.stim"

ident_1_name = paste0(i,'_',args[3])
ident_2_name = paste0(i,'_',args[4])

#<3
aaa<-table(Idents(object = combined))


if((ident_1_name %in% names(aaa)) && (ident_2_name %in% names(aaa)) ){
	if(aaa[ident_1_name]>=3 &&  aaa[ident_2_name]>=3){

	#b.interferon.response <- FindMarkers(object = combined, ident.1 = ident_1_name, ident.2 = ident_2_name, 
    #verbose = FALSE)
    b.interferon.response <- NA
    tryCatch({
    	b.interferon.response <-FindMarkers(object = combined, ident.1 = ident_1_name, ident.2 = ident_2_name ,verbose = FALSE,logfc.threshold = 0.25)
        b.interferon.response=b.interferon.response %>% filter(p_val_adj <= 0.05)
    	error = function(err){
    		cat("**********************\n",ident_1_name,ident_2_name,"\nNo features pass logfc.threshold threshold \nERROR! signed by cwt\n***********************\n")
    	}
        })
        if(length(which(! is.na(b.interferon.response) %in% 'FALSE')) == 0){
# if(! is.na(b.interferon.response)){
		head(x = b.interferon.response, n = 15)
		table_name <- paste0(args[2],"/", ident_1_name, '--vs--', ident_2_name, '_FindMarkers.xls')
		write.table(b.interferon.response, file = table_name, sep="	", col.names = NA, row.names = T, quote=F)
		i_name = gsub('/| ','_',i)
		i_name = gsub("[(.*)]","",i_name)
		features_plot_gene <- NA

		if (is.na(features_plot_gene)){
			features.plot_gene = row.names(head(b.interferon.response, 3))
		}else{
			features.plot_gene <- strsplit(features.plot_gene ,',')[[1]]
		}

		FeatureHeatmap_name = paste0(args[2],"/",i_name,'_FeatureHeatmap_name','.pdf')
		FeatureHeatmap_name_png = paste0(args[2],"/",i_name,'_FeatureHeatmap_name','.png')
		pdf(FeatureHeatmap_name)
		df<-FeaturePlot(object = combined, features = features.plot_gene, split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))
		print(df)
		dev.off()
		convert_name = paste("/usr/bin/convert  -density 300 -resize 30%",FeatureHeatmap_name,FeatureHeatmap_name_png)
		system(convert_name)



		VlnPlot_name_pdf = paste0(args[2],"/",i_name,'_VlnPlot_name','.pdf')
		VlnPlot_name_png = paste0(args[2],"/",i_name,'_VlnPlot_name','.png')
		pdf(VlnPlot_name_pdf,width=15,height=15)

		plots <- VlnPlot(object = combined, features = features.plot_gene, split.by = "stim", 
    		group.by = "celltype", pt.size = 0, combine = FALSE)

		df1 <- CombinePlots(plots = plots, ncol = 1)
		print(df1)
		dev.off()
		convert_name = paste("/usr/bin/convert  -density 300 -resize 30%",VlnPlot_name_pdf,VlnPlot_name_png)
		system(convert_name)
	}
}}
}
abc(combined)
#combined <- combined.tmp
}

#saveRDS(combined, file = "/scRNA/USER/zhouli/project/tmp/ER047.wengshijun/result/process/vis_changes_gene/Merge1/old--vs--D0/cell_Type_name_combined.rds") 


