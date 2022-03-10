library(zellkonverter)
sce <- zellkonverter::readH5AD("../../../Downloads/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad")
sce

mata <- SingleCellExperiment::colData(sce)

sce_immune <- zellkonverter::readH5AD("../../../Downloads/GTEx_8_tissues_snRNAseq_immune_atlas_071421.public_obs.h5ad")
meta <- SingleCellExperiment::colData(sce_immune)
