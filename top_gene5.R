



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

ace = run.ACTIONet(sce,k_max=strtoi(args[3]),layout_compactness=0,thread_no=30,layout_epochs = 2500)



select.top.k.features <- function(feature.enrichment.table, top.features = 3, normalize = F,
    reorder.columns = T) {

    W0 = (feature.enrichment.table)
    if (normalize == T)
        W0 = doubleNorm(W0)


    IDX = matrix(0, nrow = top.features, ncol = ncol(W0))
    VV = matrix(0, nrow = top.features, ncol = ncol(W0))
    W = (W0)
    for (i in 1:nrow(IDX)) {
        W.m = as(MWM_hungarian(W), "dgTMatrix")
  
	IDX[i, W.m@j + 1] = W.m@i + 1
        VV[i, W.m@j + 1] = W.m@x

        W[IDX[i, W.m@j + 1], ] = 0
    }

    if (reorder.columns == T) {
        feature.enrichment.table.aggregated = apply(IDX, 2, function(perm) as.numeric(Matrix::colMeans(W0[perm,
            ])))
        CC = cor(feature.enrichment.table.aggregated)
        D = as.dist(1 - CC)
        cols = seriation::get_order(seriation::seriate(D, "OLO"))
        rows = as.numeric(IDX[, cols])
    } else {
        cols = 1:ncol(W0)
        rows = unique(as.numeric(IDX))
    }
    W = feature.enrichment.table[rows, cols]

    W = t(apply(W, 1, scale))
  
    return(W)
}

plot.top.k.features <- function(feature.enrichment.table, top.features = 3, normalize = T,
    reorder.columns = T, row.title = "Archetypes", column.title = "Genes", rowPal = "black") {

    W = select.top.k.features(feature.enrichment.table, top.features = top.features, normalize = normalize)

    gradPal = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))))(100)

    M = apply(W, 1, max)

	Z = W
    #Z = Matrix::t(apply(W, 1, function(w) (w - min(w))/(max(w) - min(w))))
    # Z = t(scale(t(W))) Z[Z > 3] = 3 Z[Z < -3] = -3
    ht = Heatmap(Z, name = "Expression (scaled)", cluster_rows = F, cluster_columns = F,
        col = gradPal, row_title = row.title, column_title = column.title, column_names_gp = gpar(fontsize = 8,
            fontface = "bold"), row_names_gp = gpar(fontsize = 5, fontface = "bold",
            col = rowPal), column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        row_title_gp = gpar(fontsize = 14, fontface = "bold"), row_names_side = "left",
        rect_gp = gpar(col = "black"), row_names_max_width = unit(150, "cm"), column_names_max_height = unit(150, "cm"))

    return(ht)
}

plot.top.k.genes <- function(ace, top.genes = 5, CPal = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH",
    top.features = 3, normalize = F, reorder.columns = T, row.title = "Archetypes",
    column.title = "Genes", rowPal = "black", specificity.slot.name = "unified_feature_specificity") {
    require(ComplexHeatmap)

    feature.enrichment.table = as.matrix(rowMaps(ace)[[specificity.slot.name]])
    filtered.rows = grep(blacklist.pattern, rownames(feature.enrichment.table))
    if (length(filtered.rows) > 0)
        feature.enrichment.table = feature.enrichment.table[-filtered.rows, ]

    ht = plot.top.k.features(feature.enrichment.table, top.features = top.features,
        normalize = normalize, reorder.columns = reorder.columns, row.title = row.title,
        column.title = column.title, rowPal = rowPal)

    return(ht)

}

pdf(sprintf("Top_gene5.pdf"))
plot.top.k.genes(ace = ace, top.features = 5, normalize = F)
dev.off()



