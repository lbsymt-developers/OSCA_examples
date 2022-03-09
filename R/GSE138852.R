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
rownames(meta) <- meta$sampleID

colData(sce_138852) <- DataFrame(meta)

stopifnot(identical(rownames(meta), colnames(mat)))

sce_138852 <- scuttle::addPerFeatureQC(sce_138852)
rowData(sce_138852)
rowRanges(sce_138852)

# rowData(sce_138852)$Length <- gene.length
# rowData(sce_138852)

# library(org.Hs.eg.db)
# ensids <- rownames(sce_138852)
# cols <- c("ENSEMBL", "GENENAME")
# ens_symbol <- select(org.Hs.eg.db, keys=ensids, columns=cols, keytype="SYMBOL")
# ens_symbol <- na.omit(ens_symbol)

# NOTA, FALTA LA PARTE DE ANOTACION DE GENES Y CONVERSION A GRANGES

# sce_frankie <- sce_138852[(ens_symbol$SYMBOL),]
#
# rowData(sce_frankie)$ensids <- DataFrame(ens_symbol)
#
# # Ensembl 100.
# "http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
# library(BiocFileCache)
# path <- tempfile()
# bfc <- BiocFileCache(path, ask = FALSE)
# mm10.gtf <- bfcrpath(bfc, file.path("http://ftp.ensembl.org/pub/release-100",
#                                     "gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"))
# gene.data <- rtracklayer::import(mm10.gtf)
# gene.data <- gene.data[gene.data$type=="gene"]
# names(gene.data) <- gene.data$gene_id
# is.gene.related <- grep("gene_", colnames(mcols(gene.data)))
# mcols(gene.data) <- mcols(gene.data)[,is.gene.related]
#
# gens <- rowData(sce_frankie)$ensids$ENSEMBL
# rowRanges(sce_138852) <- gene.data[gene.data$gene_name %in% row.names(sce_138852)]
# rowRanges(sce)[1:10,]
#
# gene.data[]

#  REDUCCION DE DIMENSIONES



