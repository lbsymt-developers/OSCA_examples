library(BiocFileCache)
path <- tempfile()
bfc <- BiocFileCache(path, ask = FALSE)
url <- file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
                 "GSE138nnn/GSE138852/suppl",
                 "GSE138852%5Fcounts%2Ecsv%2Egz")

grubman.fname <- bfcrpath(bfc, url)
local.name <- URLdecode(basename(url))
unlink(local.name)
if (.Platform$OS.type=="windows") {
  file.copy(grubman.fname, local.name)
} else {
  file.symlink(grubman.fname, local.name)
}

aa <- data.table::fread("GSE138852_counts.csv.gz")
colnames(aa)[1] <- "gene"
rownames(aa) <- aa$gene
genes <- aa$gene
aa <- aa[,-1]
mat <- as.matrix(aa)
row.names(mat) <- genes

library(Matrix)
sparse.mat <- as(mat, "sparseMatrix")
dim(sparse.mat)

saveRDS(sparse.mat, "data/GSE138852_sparseMatrix.rds")



# sparse.mat <- readSparseCounts(data.table::fread("GSE138852_counts.csv.gz"))
# Este no se pudo, hay que indagar la causa


# Creamos el objeto singlecellexperiment
library(SingleCellExperiment)
mat <- readRDS("data/GSE138852_sparseMatrix.rds")
sce_138852 <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = mat))

assay(sce_138852, "counts")
sce_138852 <- scuttle::logNormCounts(sce_138852)
sce_138852
assay(sce_138852, "logcounts")

meta <- readr::read_delim("data/scRNA_metadata.tsv")


