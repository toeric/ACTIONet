



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


print(sprintf("ACTION good !"))

# clustering
clustering.resolution <- as.numeric(args[4])
print(sprintf("resolution: %f", clustering.resolution))

clusters = Leiden.clustering(ACTIONet_results, resolution_parameter = clustering.resolution)

print(sprintf("ACTION good2 !"))


pdf(sprintf("cluter_result.pdf"))
plot.ACTIONet(ACTIONet_results, clusters)
dev.off()

print(sprintf("ACTION good3 !"))

# correct cell annotations
#min.cell.fraction = 0.0005
#ACTIONet.out = correct.cell.annotations(clusters, min.cell.fraction = min.cell.fraction, adjust.levels=TRUE)

ACTIONet_results = compute.cluster.feature.specificity(ACTIONet_results, clusters, output.slot.name = "Leiden_specificity")

ACTIONet_results


spe = ACTIONet_results$Leiden_specificity_feature_specificity
spe

write.table(as.matrix(spe), file="specificity.txt", quote=FALSE,sep=",",col.names=TRUE)

