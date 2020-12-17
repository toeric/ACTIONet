args<-commandArgs(TRUE)
require(ACTIONet)
require(SingleCellExperiment)
library(dplyr)
#library(Seurat)
require(ComplexHeatmap)

COLOR_PHASE = c("#FC6556","#29C744","#6169FF","#F9F330","#FFA832","#FF9FC0","#AB7DFF","#BFC9C3")

setwd(args[1])
part_name <- strsplit(args[2],"[-.]") %>% unlist()
part_name <- part_name[4]


#----------------------------------------------#
# Read the reduced SCE object and run ACTIONet #
#----------------------------------------------#
sce = readRDS(args[2])

ACTIONet_results = run.ACTIONet(sce,k_max=strtoi(args[3]),layout_compactness=0,thread_no=30,layout_epochs = 2500)

print(sprintf("ACTION good !"))


gene_exp = as.matrix(rowMaps(ACTIONet_results)[["unified_feature_profile"]])
write.table(gene_exp, file=sprintf("gene_exp.txt"))
