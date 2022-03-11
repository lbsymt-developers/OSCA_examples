# Lo demostramos utilizando un pequeño conjunto de datos de scRNA-seq de Lun et al. (2017),
# que se proporciona sin control de calidad previo para que podamos aplicar nuestros propios procedimientos.

#--- loading ---#
library(scRNAseq)
sce.416b <- LunSpikeInData(which="416b")
sce.416b$block <- factor(sce.416b$block)
sce.416b

# Identifying the mitochondrial transcripts in our SingleCellExperiment.
location <- rowRanges(sce.416b)
is.mito <- any(seqnames(location)=="MT")


# Para cada célula, calculamos estas métricas de control de calidad utilizando la función perCellQCMetrics()
# del paquete scater (McCarthy et al. 2017). La columna de la suma contiene el recuento total para cada célula y
# la columna detectada contiene el número de genes detectados.
# La columna subsets_Mito_percent contiene el porcentaje de lecturas mapeadas a transcripciones mitocondriales.
# Por último, la columna altexps_ERCC_percent contiene el porcentaje de lecturas asignadas a transcripciones ERCC.

library(scuttle)
df <- perCellQCMetrics(sce.416b, subsets=list(Mito=is.mito))
summary(df$sum)
summary(df$detected)
summary(df$subsets_Mito_percent)
summary(df$altexps_ERCC_percent)


# Como alternativa, los usuarios pueden preferir utilizar la función addPerCellQC().
# Esta función calcula y añade las estadísticas QC por celda a los colData del objeto SingleCellExperiment,
# lo que nos permite conservar toda la información relevante en un único objeto para su posterior manipulación.

sce.416b <- addPerCellQCMetrics(sce.416b, subsets=list(Mito=is.mito))
colnames(colData(sce.416b))


# El enfoque más sencillo para identificar las células de baja calidad consiste en aplicar umbrales fijos
# a las métricas de control de calidad. Por ejemplo, podríamos considerar que las células son de baja calidad
# si tienen un tamaño de biblioteca inferior a 100.000 lecturas; expresan menos de 5.000 genes;
# tienen proporciones de espigas superiores al 10%; o tienen proporciones mitocondriales superiores al 10%.

qc.lib <- df$sum < 1e5
qc.nexprs <- df$detected < 5e3
qc.spike <- df$altexps_ERCC_percent > 10
qc.mito <- df$subsets_Mito_percent > 10
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito

# Summarize the number of cells removed for each reason.
DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
          SpikeProp=sum(qc.spike), MitoProp=sum(qc.mito), Total=sum(discard))




