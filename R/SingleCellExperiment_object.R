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

# Con esto añadimos las cuentas manualmente con modificaciones
counts_100 <- counts(sce) + 100
assay(sce, "counts_100") <- counts_100 # assign a new entry to assays slot
assays(sce) # new assay has now been added.


# Para obtener los ensayos en el objeto singleCellExperiment se usa:
assays(sce)
# Only keeping the first two assays
assays(sce) <- assays(sce)[1:2]
sce

assayNames(sce)
names(assays(sce)) # same result, but slightly less efficient.


# Añadiendo metadata
# Downloading the SDRF file containing the metadata.
lun.sdrf <- bfcrpath(bfc, file.path("https://www.ebi.ac.uk/arrayexpress/files",
                                    "E-MTAB-5522/E-MTAB-5522.sdrf.txt"))

coldata <- read.delim(lun.sdrf, check.names=FALSE)

# Only keeping the cells involved in the count matrix in 'mat'.
coldata <- coldata[coldata[,"Derived Array Data File"]=="counts_Calero_20160113.tsv",]

# Only keeping interesting columns, and setting the library names as the row names.
coldata <- DataFrame(
  genotype=coldata[,"Characteristics[genotype]"],
  phenotype=coldata[,"Characteristics[phenotype]"],
  spike_in=coldata[,"Factor Value[spike-in addition]"],
  row.names=coldata[,"Source Name"]
)

coldata

sce <- SingleCellExperiment(assays = list(counts=mat), colData=coldata)
sce

# Formas de añadir Metadata
colData(sce)
head(sce$Factor.Value.phenotype.)

sce <- SingleCellExperiment(list(counts=mat))
colData(sce) <- coldata

sce <- SingleCellExperiment(list(counts=mat))
sce$phenotype <- coldata$phenotype
colData(sce)
sce

# Para estar seguros de la consistencia de los nombres
stopifnot(identical(rownames(coldata), colnames(mat)))

# Esto nos ayuda a agregar una metrica de calidad como metadata
sce <- scuttle::addPerCellQC(sce)
colData(sce)

# Almacenamos la anotación a nivel de característica en la ranura rowData,
# un DataFrame donde cada fila corresponde a un gen y contiene anotaciones como la longitud de la transcripción o el símbolo del gen.
# Podemos obtener y establecer esto manualmente con la función rowData(), como se muestra a continuación:
rowData(sce)$Length <- gene.length
rowData(sce)

sce <- scuttle::addPerFeatureQC(sce)
rowData(sce)

# Furthermore, there is a special rowRanges slot to hold genomic coordinates in the form of a GRanges or GRangesList.
rowRanges(sce) # empty

# La forma en que rellenamos los rowRanges depende del organismo y de la anotación utilizada durante el alineamiento y la cuantificación.
# En este caso, tenemos identificadores de Ensembl, por lo que podríamos utilizar rtracklayer
# para cargar un GRanges desde un archivo GTF que contenga la anotación de Ensembl utilizada en este conjunto de datos:
gene.data <- rtracklayer::import(mm10.gtf)

# Cleaning up the object.
gene.data <- gene.data[gene.data$type=="gene"]
names(gene.data) <- gene.data$gene_id
is.gene.related <- grep("gene_", colnames(mcols(gene.data)))
mcols(gene.data) <- mcols(gene.data)[,is.gene.related]

rowRanges(sce) <- gene.data[rownames(sce)]
rowRanges(sce)[1:10,]

# OTROS METADATOS



