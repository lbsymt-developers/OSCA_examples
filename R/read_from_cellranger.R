library(DropletUtils)
sce <- read10xCounts("tenx-2.1.0-pbmc4k/filtered_gene_bc_matrices/GRCh38")
sce

dirA <- "tenx-2.1.0-pbmc4k/filtered_gene_bc_matrices/GRCh38"
dirB <- "tenx-2.1.0-pbmc4k-2/filtered_gene_bc_matrices/GRCh38"
sce <- read10xCounts(c(dirA, dirB))
sce
