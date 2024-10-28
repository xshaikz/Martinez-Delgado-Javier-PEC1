library(readxl)
library(SummarizedExperiment)

# Creamos un archivo temporal en el que guardar los datos y los cargamos en un dataframe para trabajar en ellos
url = "https://github.com/nutrimetabolomics/metaboData/raw/main/Datasets/2018-Phosphoproteomics/TIO2%2BPTYR-human-MSS%2BMSIvsPD.XLSX"
temporal = tempfile(fileext = ".xlsx")
download.file(url, destfile = temporal, mode = "wb")
data = read_excel(temporal)

# Comprobamos que los datos se han cargado bien
head(data)

# Creamos la matriz con los valores de expresión
matriz_expresion = as.matrix(data[,5:16])

columnas_expresion = colnames(matriz_expresion)

# Sacamos como metadatos si es MSS o PD
vector_vacio = character(length(columnas_expresion))
for (i in seq_along(columnas_expresion)) {
  nombre_columna = columnas_expresion[i]
  parte = strsplit(nombre_columna, "_")[[1]]
  condicion = parte[length(parte)] 
  vector_vacio[i] = condicion
}

# Juntamos el nombre de las columnas con los datos de si es MSS o PD para crear los metadatos de la muestra
metadatos_muestra = data.frame(
  id_muestra = columnas_expresion,
  condicion = vector_vacio
)

# Seleccionamos las demás variables que no son de la muestra para crear los metadatos de las variables
metadatos_variables = data[,c("SequenceModifications","Accession","Description","Score","CLASS","PHOSPHO")]

# Creamos el objeto SummarizedExperiment
se = SummarizedExperiment(
  assays = list(counts = matriz_expresion),
  rowData = metadatos_variables,
  colData = metadatos_muestra
)

# Resumen estadístico
summary(assay(se, "counts"))

# Histogramas
library(ggplot2)
log_expr_data = log2(expr_data + 1)
log_expr_data_long = reshape2::melt(log_expr_data, variable.name = "Muestra", value.name = "Intensidad_Log")

ggplot(log_expr_data_long, aes(x = Intensidad_Log, color = Muestra)) +
  geom_density() +
  labs(title = "Distribución de intensidad de expresión por muestra", 
       x = "Intensidad (Log2)", y = "Densidad") +
  theme_minimal()


boxplot(log2(assay(se, "counts")), main="Boxplot de intensidades de expresión por muestra", 
        xlab="Muestras", ylab="Intensidad", las=1, col="lightblue")

# Definimos la tolerancia para detectar una varianza baja
tolerancia = 1e-8
expr_matrix_filtered_rows = expr_matrix[apply(expr_matrix, 1, var) > tolerancia, ]
expr_matrix_filtered = expr_matrix_filtered_rows[, apply(expr_matrix_filtered_rows, 2, var) > tolerancia]

# PCA
pca = prcomp(t(expr_matrix_filtered), scale. = TRUE)
pca_df = as.data.frame(pca$x)
pca_df$condition = colData(se)$condition[colnames(expr_matrix_filtered)]

ggplot(pca_df, aes(x = PC1, y = PC2, color = condicion)) +
  geom_point(size = 3) +
  labs(title = "PCA de las muestras (filtrado de filas y columnas)", x = "PC1", y = "PC2") +
  theme_minimal()

save(se, file = "objeto_contenedor.rda")

expression_data = as.data.frame(assay(se, "counts"))
write.csv(expression_data, "datos_expresion.csv", row.names = TRUE)


library(knitr)
sample_metadata = as.data.frame(colData(se))
sample_md = knitr::kable(sample_metadata, format = "markdown")
cat(sample_md, file = "metadatos_muestra.md", sep = "\n")

variable_metadata = as.data.frame(rowData(se))
variable_md = knitr::kable(variable_metadata, format = "markdown")
cat(variable_md, file = "metadatos_variables.md", sep = "\n")
