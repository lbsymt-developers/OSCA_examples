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

sce_138852 <- scater::logNormCounts(sce_138852)
sce_138852 <- scater::runPCA(sce_138852)
dim(reducedDim(sce_138852, "PCA"))

sce_138852 <- scater::runTSNE(sce_138852, perplexity = 0.1)
head(reducedDim(sce_138852, "TSNE"))

sce_138852 <- scater::runUMAP(sce_138852)
# reducedDim(sce, "UMAP_uwot") <- u
# reducedDims(sce) # Now stored in the object.
# head(reducedDim(sce, "UMAP_uwot"))

save(sce_138852, file = "data/sce_138852.RData")

load("data/sce_138852.RData")

#########################################################
#########################################################
#########################################################
# CONTROL DE CALIDAD
library(scuttle)
location <- rowRanges(sce_138852)
is.mito <- any(seqnames(location)=="MT")
df <- perCellQCMetrics(sce_138852, subsets=list(Mito=is.mito))
summary(df$sum)
summary(df$detected)
summary(df$subsets_Mito_percent)
summary(df$altexps_ERCC_percent)

sce_138852 <- addPerCellQCMetrics(sce_138852, subsets=list(Mito=is.mito))
colnames(colData(sce_138852))

qc.lib <- df$sum < 1e5
qc.nexprs <- df$detected < 5e3
qc.spike <- df$altexps_ERCC_percent > 10
qc.mito <- df$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

reasons <- perCellQCFilters(df,
                            sub.fields=c("subsets_Mito_percent"))
colSums(as.matrix(reasons))

DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
          SpikeProp=sum(qc.spike), MitoProp=sum(qc.mito), Total=sum(discard))

stats <- cbind(log10(df$sum), log10(df$detected),
               df$subsets_Mito_percent, df$altexps_ERCC_percent)

colData(sce_138852) <- cbind(colData(sce_138852), df)
sce_138852$batch <- factor(sce_138852$batch)
sce_138852$batchCond <- ifelse(grepl("induced", sce_138852$phenotype),
                             "induced", "wild type")

sce_138852$discard <- reasons$discard

library(scater)
tiff("images/diagnostic_plot.tiff", height = 30, width = 25, units='cm',
     compression = "lzw", res = 300)
gridExtra::grid.arrange(
  plotColData(sce_138852, x="batch", y="sum", colour_by="discard",
              other_fields="batchCond") + facet_wrap(~batchCond) +
    scale_y_log10() + ggtitle("Total count") +
    theme(text = element_text(size = 20)),
  plotColData(sce_138852, x="batch", y="detected", colour_by="discard",
              other_fields="batchCond") + facet_wrap(~batchCond) +
    scale_y_log10() + ggtitle("Detected features") +
    theme(text = element_text(size = 20)),
  plotColData(sce_138852, x="batch", y="subsets_Mito_percent",
              colour_by="discard", other_fields="batchCond") +
    facet_wrap(~batchCond) + ggtitle("Mito percent") +
    theme(text = element_text(size = 20)),
  ncol=1
)
dev.off()

sce_138852 <- addPerCellQC(sce_138852,
                           subsets=list(Mt=rowData(sce_138852)$featureType=="mito"))
qc <- quickPerCellQC(colData(sce_138852),
                     sub.fields=c("subsets_Mt_percent"))
sce_138852$discard <- qc$discard
plotColData(sce_138852, x="sum", y="subsets_Mt_percent", colour_by="discard")

tiff("images/UMAP_gse138852.tiff", height = 30, width = 25, units='cm',
     compression = "lzw", res = 300)
plotUMAP(sce_138852, colour_by="cellType", shape_by = "batchCond") +
  theme(text = element_text(size = 20))
dev.off()

library(scran)
g <- buildSNNGraph(sce_138852, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
colLabels(sce_138852) <- factor(clust)
plotUMAP(sce_138852, colour_by="label")
markers <- scran::findMarkers(sce_138852, pval.type="some", direction="up")
table(colLabels(sce_138852))
marker.set <- markers[["11"]]
as.data.frame(marker.set[1:30,1:3])

tiff("images/expression_genes.tiff", height = 30, width = 25, units='cm',
     compression = "lzw", res = 300)
plotExpression(sce_138852, features=c("APOE", "TOMM40",
                                    "APOC1", "CLU"), x="label", colour_by="label")+
  theme(text = element_text(size = 20))
dev.off()


tiff("images/UMAP_labels.tiff", height = 30, width = 25, units='cm',
     compression = "lzw", res = 300)
plotUMAP(sce_138852, colour_by="label")+
  theme(text = element_text(size = 20))
dev.off()

# VARIANZA

set.seed(1001)
dec.pbmc <- modelGeneVarByPoisson(sce_138852)
top.pbmc <- getTopHVGs(dec.pbmc, prop=0.1)
plot(dec.pbmc$mean, dec.pbmc$total, pch=16, cex=0.5,
     xlab="Mean of log-expression", ylab="Variance of log-expression")+
  theme(text = element_text(size = 20))
curfit <- metadata(dec.pbmc)
curve(curfit$trend(x), col='dodgerblue', add=TRUE, lwd=2)
