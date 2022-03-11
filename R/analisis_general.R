# Aquí utilizamos el conjunto de datos de retina basado en gotas de Macosko et al. (2015),
# proporcionado en el paquete scRNAseq. Se parte de una matriz de recuento y se termina con clusters en preparación para la interpretación biológica.
# Hay flujos de trabajo similares en forma abreviada en partes posteriores del libro.


library(scRNAseq)
sce <- MacoskoRetinaData()

# Quality control (using mitochondrial genes).
library(scater)
is.mito <- grepl("^MT-", rownames(sce))
qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))
filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")
sce <- sce[, !filtered$discard]

# Normalization.
sce <- logNormCounts(sce)

# Feature selection.
library(scran)
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, prop=0.1)

# PCA.
library(scater)
set.seed(1234)
sce <- runPCA(sce, ncomponents=25, subset_row=hvg)

# Clustering.
library(bluster)
colLabels(sce) <- clusterCells(sce, use.dimred='PCA',
                               BLUSPARAM=NNGraphParam(cluster.fun="louvain"))

# Visualization.
sce <- runUMAP(sce, dimred = 'PCA')
plotUMAP(sce, colour_by="label")

# Marker detection.
markers <- findMarkers(sce, test.type="wilcox", direction="up", lfc=1)


# Aquí utilizamos el conjunto de datos Smart-seq2 del páncreas de Segerstolpe et al. (2016),
# de nuevo proporcionado en el paquete scRNAseq. Se parte de una matriz de recuento y se termina con clusters
# con algunos ajustes adicionales para eliminar los efectos de lote no interesantes entre individuos.
# Tenga en cuenta que un análisis más elaborado del mismo conjunto de datos con justificaciones
# para cada paso está disponible en el capítulo 8 del flujo de trabajo.


sce <- SegerstolpePancreasData()

# Quality control (using ERCCs).
qcstats <- perCellQCMetrics(sce)
filtered <- quickPerCellQC(qcstats, percent_subsets="altexps_ERCC_percent")
sce <- sce[, !filtered$discard]

# Normalization.
sce <- logNormCounts(sce)

# Feature selection, blocking on the individual of origin.
dec <- modelGeneVar(sce, block=sce$individual)
hvg <- getTopHVGs(dec, prop=0.1)

# Batch correction.
library(batchelor)
set.seed(1234)
sce <- correctExperiments(sce, batch=sce$individual,
                          subset.row=hvg, correct.all=TRUE)

# Clustering.
colLabels(sce) <- clusterCells(sce, use.dimred='corrected')

# Visualization.
sce <- runUMAP(sce, dimred = 'corrected')
gridExtra::grid.arrange(
  plotUMAP(sce, colour_by="label"),
  plotUMAP(sce, colour_by="individual"),
  ncol=2
)

# Marker detection, blocking on the individual of origin.
markers <- findMarkers(sce, test.type="wilcox", direction="up", lfc=1)
