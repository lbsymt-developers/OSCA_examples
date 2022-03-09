library(zellkonverter)
demo <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")
sce <- readH5AD(demo)
sce

library(LoomExperiment)
demo <- system.file("extdata", "L1_DRG_20_example.loom", package = "LoomExperiment")
scle <- import(demo, type="SingleCellLoomExperiment")
scle
