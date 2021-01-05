



args<-commandArgs(TRUE)
require(ACTIONet)
require(SingleCellExperiment)
library(dplyr)
#library(Seurat)

COLOR_PHASE = c("#FC6556","#29C744","#6169FF","#F9F330","#FFA832","#FF9FC0","#AB7DFF","#BFC9C3")

setwd(args[1])
part_name <- strsplit(args[2],"[-.]") %>% unlist()
part_name <- part_name[4]


#----------------------------------------------#
# Read the reduced SCE object and run ACTIONet #
#----------------------------------------------#
sce = readRDS(args[2])

ACTIONet_results = run.ACTIONet(sce,k_max=strtoi(args[3]),layout_compactness=0,thread_no=30,layout_epochs = 2500)

ACTIONet_results

str(ACTIONet_results)

#gene_dif = as.matrix(rowMaps(ACTIONet_results)[["unified_feature_specificity"]])
#write.table(gene_dif, file=sprintf("gene_dif.txt"))


#gene_exp = as.matrix(rowMaps(ACTIONet_results)[["unified_feature_profile"]])
#write.table(gene_exp, file=sprintf("gene_exp.txt"))

print(sprintf("ACTION done !"))

# clustering
clustering.resolution <- as.numeric(args[4])
print(sprintf("resolution: %f", clustering.resolution))
clusters = Leiden.clustering(ACTIONet_results, resolution_parameter = clustering.resolution)

print(sprintf("Cluster done !"))

pdf(sprintf("cluter_result.pdf"))
plot.ACTIONet(ACTIONet_results, clusters)
dev.off()


ACTIONet_results = compute.cluster.feature.specificity(ACTIONet_results, clusters, output.slot.name = "Leiden_specificity")


print(sprintf("specificity done !"))


spe = ACTIONet_results$Leiden_specificity_feature_specificity

write.table(as.matrix(spe), file="specificity.txt", quote=FALSE,sep=",",col.names=TRUE)



