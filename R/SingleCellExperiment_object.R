library(BiocFileCache)
bfc <- BiocFileCache("raw_data", ask = FALSE)
lun.zip <- bfcrpath(bfc, file.path("https://www.ebi.ac.uk/arrayexpress/files",
                                   "E-MTAB-5522/E-MTAB-5522.processed.1.zip"))
unzip(lun.zip, exdir=tempdir())

mat <- read.delim(file.path(tempdir(), "counts_Calero_20160113.tsv"),
                  header=TRUE, row.names=1, check.names=FALSE)

# Only considering endogenous genes for now.
spike.mat <- mat[grepl("^ERCC-", rownames(mat)),]
mat <- mat[grepl("^ENSMUSG", rownames(mat)),]

# Splitting off the gene length column.
gene.length <- mat[,1]
mat <- as.matrix(mat[,-1])

dim(mat)

library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = mat))
sce

assay(sce, "counts")
counts(sce)
mat2 <- counts(sce)
sce <- scuttle::logNormCounts(sce)
sce
