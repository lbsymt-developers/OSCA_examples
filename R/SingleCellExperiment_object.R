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

"http://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/"
# The paper uses Ensembl 82.
library(BiocFileCache)
mm10.gtf <- bfcrpath(bfc, file.path("http://ftp.ensembl.org/pub/release-82",
                                    "gtf/mus_musculus/Mus_musculus.GRCm38.82.gtf.gz"))
gene.data <- rtracklayer::import(mm10.gtf)

# Cleaning up the object.
gene.data <- gene.data[gene.data$type=="gene"]
names(gene.data) <- gene.data$gene_id
is.gene.related <- grep("gene_", colnames(mcols(gene.data)))
mcols(gene.data) <- mcols(gene.data)[,is.gene.related]

rowRanges(sce) <- gene.data[rownames(sce)]
rowRanges(sce)[1:10,]

# OTROS METADATOS

# Para agregar meytadatos de genes que nos interesen
my_genes <- c("gene_1", "gene_5")
metadata(sce) <- list(favorite_genes = my_genes)
metadata(sce)

# Y se puede agregar mas información con el simbolo $
your_genes <- c("gene_4", "gene_8")
metadata(sce)$your_genes <- your_genes
metadata(sce)


# SUBSETTING Y COMBINACIÓN

first.10 <- sce[,1:10]
ncol(counts(first.10)) # only 10 columns.

# Del mismo modo, si sólo quisiéramos celdas de tipo salvaje,
# podríamos subconjuntar nuestro objeto sce basándonos en sus entradas de colData():
wt.only <- sce[, sce$phenotype == "wild type phenotype"]
ncol(counts(wt.only))
colData(wt.only)

# La misma lógica se aplica a rowData().
# Digamos que sólo queremos conservar los genes codificadores de proteínas:
coding.only <- sce[rowData(sce)$gene_biotype == "protein_coding",]
nrow(counts(coding.only))
rowData(coding.only)

# Por el contrario, si tuviéramos que combinar varios objetos SingleCellExperiment,
# la clase se encargaría de combinar tanto los valores de expresión como la anotación asociada de forma coherente.
# Podemos utilizar cbind() para combinar objetos por columna,
# suponiendo que todos los objetos implicados tienen los mismos valores de anotación de fila y campos de anotación de columna compatibles.

sce2 <- cbind(sce, sce)
ncol(counts(sce2)) # twice as many columns
colData(sce2) # twice as many rows

# Del mismo modo, podemos utilizar rbind() para combinar objetos por fila,
# suponiendo que todos los objetos tienen los mismos valores de anotación de columna y campos de anotación de fila compatibles.
sce2 <- rbind(sce, sce)
nrow(counts(sce2)) # twice as many rows
rowData(sce2) # twice as many rows


# SINGLE-CELL-SPECIFIC-FIELDS

# Reduccion de dimensiones
sce <- scater::logNormCounts(sce)
sce <- scater::runPCA(sce)
dim(reducedDim(sce, "PCA"))

sce <- scater::runTSNE(sce, perplexity = 0.1)
head(reducedDim(sce, "TSNE"))

# Podemos ver los nombres de todas nuestras entradas en la ranura reducedDims a través del accesorio,
# reducedDims(). Tenga en cuenta que esto es plural y devuelve una lista de todos los resultados, mientras que reducedDim() sólo devuelve un único resultado.

reducedDims(sce)

# Para agregar el clustering o reduccion manualmente se utiliza lo siguiente:
u <- uwot::umap(t(logcounts(sce)), n_neighbors = 2)
reducedDim(sce, "UMAP_uwot") <- u
reducedDims(sce) # Now stored in the object.
head(reducedDim(sce, "UMAP_uwot"))


# Factores de tamaño
sce <- scran::computeSumFactors(sce)
summary(sizeFactors(sce))

sizeFactors(sce) <- scater::librarySizeFactors(sce)
summary(sizeFactors(sce))

# La función colLabels() nos permite obtener o establecer un vector o factor de etiquetas por célula,
# que suelen corresponder a las agrupaciones asignadas por el clustering no supervisado (Capítulo Básico 5)
# o a las identidades de tipo de célula predichas por los algoritmos de clasificación (Capítulo Básico 7).

colLabels(sce) <- scran::clusterCells(sce, use.dimred="PCA")
table(colLabels(sce))

