library(scRNAseq)
sce.zeisel <- ZeiselBrainData()

library(scater)
sce.zeisel <- aggregateAcrossFeatures(sce.zeisel,
                                      id=sub("_loc[0-9]+$", "", rownames(sce.zeisel)))
# library(org.Mm.eg.db)
# rowData(sce.zeisel)$Ensembl <- mapIds(org.Mm.eg.db,
#                                       keys=rownames(sce.zeisel), keytype="SYMBOL", column="ENSEMBL")




