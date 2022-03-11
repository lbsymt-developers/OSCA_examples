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


# En este caso, suponemos que la mayor parte del conjunto de datos está formada por celdas de alta calidad.
# A continuación, identificamos las celdas que son valores atípicos para las distintas métricas de control de calidad,
# basándonos en la desviación absoluta media (MAD) del valor medio de cada métrica en todas las celdas. Por defecto,
# consideramos que un valor es un valor atípico si se aleja más de 3 DAM de la mediana en la dirección "problemática".
# Esto está motivado por el hecho de que este filtro retendrá el 99% de los valores no atípicos que siguen una distribución normal.
# Lo demostramos utilizando la función perCellQCFilters() en las métricas de control de calidad del conjunto de datos 416B.

reasons <- perCellQCFilters(df,
                            sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
colSums(as.matrix(reasons))

# Esta función identificará las celdas con tamaños de bibliotecas transformados logarítmicamente que estén más de 3 DAM por debajo de la mediana.
# Se utiliza una transformación logarítmica para mejorar la resolución en valores pequeños cuando type="lower" y para evitar umbrales negativos que no tendrían sentido para una métrica no negativa.
# Además, no es infrecuente que la distribución de los tamaños de las bibliotecas muestre una cola derecha pesada;
# la transformación logarítmica evita la inflación de la MAD de una manera que podría comprometer la detección de valores atípicos en la cola izquierda.
# (En general, hace que la distribución parezca más normal para justificar el razonamiento del 99% mencionado anteriormente).
# La función también hará lo mismo para el número de genes expresados transformado en logaritmo.

# perCellQCFilters() también identificará valores atípicos para las métricas basadas en la proporción especificadas en los argumentos sub.fields=.
# Estas distribuciones a menudo muestran una cola derecha pesada, pero a diferencia de las dos métricas anteriores,
# es la propia cola derecha la que contiene las supuestas células de baja calidad.
# Por lo tanto, no realizamos ninguna transformación para reducir la cola, sino que esperamos que las celdas de la cola se identifiquen como grandes valores atípicos.
# (Aunque teóricamente es posible obtener un umbral sin sentido por encima del 100%, esto es lo suficientemente raro como para no ser una preocupación práctica).

# Una celda que es un valor atípico para cualquiera de estas métricas se considera de baja calidad y se descarta.
# Esto se recoge en la columna de descartes, que puede utilizarse para un filtrado posterior

summary(reasons$discard)
attr(reasons$low_lib_size, "thresholds")
attr(reasons$low_n_features, "thresholds")

# Con esta estrategia, los umbrales se adaptan tanto a la ubicación como a la dispersión de la distribución de valores para una métrica determinada.
# Esto permite que el procedimiento de control de calidad se ajuste a los cambios en la profundidad de la secuenciación,
# la eficiencia de la captura de ADNc, el contenido mitocondrial, etc., sin requerir ninguna intervención del usuario o experiencia previa.


# OTRAS APROXIMACIONES

# Otra estrategia consiste en identificar los valores atípicos en un espacio de alta dimensión basándose en las métricas de control de calidad de cada celda.
# Utilizamos métodos de robustbase para cuantificar la "perificidad" de cada celda en función de sus métricas de control de calidad y, a continuación,
# utilizamos isOutlier() para identificar las celdas de baja calidad que presentan niveles inusualmente altos de perificidad.

stats <- cbind(log10(df$sum), log10(df$detected),
               df$subsets_Mito_percent, df$altexps_ERCC_percent)

library(robustbase)
outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
multi.outlier <- isOutlier(outlying, type = "higher")
summary(multi.outlier)

# Para completar, observamos que los valores atípicos también pueden identificarse a partir de los perfiles de expresión génica,
# en lugar de las métricas de control de calidad.
# Consideramos que se trata de una estrategia arriesgada, ya que puede eliminar células de alta calidad en poblaciones poco frecuentes.


# PLOTS DE DIAGNOSTICO


