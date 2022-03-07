library(BiocFileCache)
path <- tempfile()
bfc <- BiocFileCache(path, ask = FALSE)
bfc <- BiocFileCache(ask=FALSE)
url <- file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
                 "GSE85nnn/GSE85241/suppl",
                 "GSE85241%5Fcellsystems%5Fdataset%5F4donors%5Fupdated%2Ecsv%2Egz")
# El link de arriba viene desde GEO Omnibus, al poner el cursor en `(ftp)`

# Hacer un enlace simbólico para que el código posterior pueda
# fingir que hemos descargado el archivo en el directorio local.
muraro.fname <- bfcrpath(bfc, url) #Esta clase representa la ubicación de los archivos almacenados en el disco. Utilice el valor de retorno para añadir y recuperar archivos que persisten a través de las sesiones.
local.name <- URLdecode(basename(url)) #Funciones para codificar o descodificar porcentualmente los caracteres de las URL.
unlink(local.name)
if (.Platform$OS.type=="windows") {
  file.copy(muraro.fname, local.name)
} else {
  file.symlink(muraro.fname, local.name)
}

# Descargar como matriz
mat <- as.matrix(read.delim("GSE85241_cellsystems_dataset_4donors_updated.csv.gz"))
dim(mat) # number of rows, number of columns


# Descargar como matriz esparcida
library(scuttle)
sparse.mat <- readSparseCounts("GSE85241_cellsystems_dataset_4donors_updated.csv.gz")
dim(sparse.mat)

# We can see that it uses less memory compared to 'mat'.
object.size(sparse.mat)
object.size(mat)
