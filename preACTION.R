args<-commandArgs(TRUE)
require(ACTIONet)
require(Seurat)
require(SingleCellExperiment)
require(scater)

setwd(dirname(args[1]))
print(sprintf("output folder: %s", dirname(args[1])))


sc = readRDS(file.path(dirname(args[1]), args[1]))
sc

sce <- as.SingleCellExperiment(sc)
sce

batch_attr = interaction(sce$source, sce$Phase)

table(batch_attr)


corrected.sce = reduce.and.batch.correct.ace.Harmony(sce, batch_attr = batch_attr, reduced_dim = 50)

saveRDS(corrected.sce, file=sprintf("sce_%s_reducedim50_batchcorrected.rds",args[2]))
