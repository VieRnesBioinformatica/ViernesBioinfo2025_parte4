#How to PCA 
#Script to inspect RNA-seq data counts 
#https://github.com/paulinapglz99

#I used the reference to build a PCA https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_03_Step_By_Step_PCA.pdf

#Instalación de paquetes necesarios --- ---
#Solo debe ejecutarse la primera vez si no tienes los paquetes instalados.

install.packages("BiocManager")
install.packages("pacman")

#Instalar los paquetes GEOquery y Biobase
BiocManager::install(c("GEOquery", "Biobase"), force = TRUE)

pacman::p_load("dplyr",         #Para manejar los datos
               "ggplot2",       #Para graficar
               "stringr",       #Para manejar los nombres de genes
               "GEOquery",      #Para descargar los datos
               "Biobase",       #Para descargar los datos
               "vroom",         #Para descargar los datos
               "tibble")       

#Descarga del conjunto de datos de expresion de GEO (GSE119834) (Pediatric Glioblastoma)

expresion <- vroom("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE119834&format=file&file=GSE119834_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
dim(expresion)

#Descarga el objeto GEO para obtener la metadata
gse <- getGEO("GSE119834")

metadata <- pData(gse[[1]])
dim(metadata)

#En el caso de GSE119834, los controles son 'Neural Stem Cells (NSCs)'.

table(metadata$source_name_ch1)

#Filtrar metadata

metadata <- metadata[metadata$geo_accession %in% colnames(expresion), ]
dim(metadata)

#Prepare data

matrix <- expresion %>% 
  column_to_rownames("GeneID") %>% 
  as.matrix() %>% 
  t()  #transpose the matrix so that rows = samples and columns = variables

#Look at the first 10 rows and first 5 columns of the matrix
matrix[1:10, 1:10]

#Ahora la Magia!

#PCA prcomp
#prcomp() asume que cada fila es una observación (sample) y cada columna una variable (feature).

pca <- prcomp(matrix,
              retx = TRUE, #valores de los puntos proyectados a la recta (los “scores”). Se almacenan en pca$x.
              center = TRUE, #resta la media de cada variable (centra)
              scale. = FALSE) #dividir por la varianza unitaria. FALSE: se usa cuando todas las variables están en la misma escala o unidades comparables (por ejemplo, conteos normalizados).

#Al ejecutar prcomp(), obtienes un objeto de clase "prcomp" con varios elementos:

names(pca)
# [1] "sdev" "rotation" "center" "scale" "x"

#Las desviaciones estándar de cada componente principal.
#Sirven para calcular la varianza explicada por cada PC. 
pca$sdev

#La matriz de carga o loadings: los coeficientes que indican cuánto contribuye cada variable original a cada componente principal.
#Cada columna es un componente principal; cada fila, una variable original.
pca$rotation[1:10, 1:10]

#Los valores medios de cada variable que fueron restados cuando center = TRUE.
pca$center

#Los factores de escala aplicados si scale. = TRUE.
#(Si no escalaste, contiene NULL.)

pca$scale

#Las coordenadas de tus observaciones proyectadas en el nuevo espacio de componentes principales.

pca$x[1:10, 1:10]

#PCA to table --- --- 

pca_df <- pca$x %>%
  as.data.frame() %>%
  rownames_to_column(var = 'IndividualID')

pca_df[1:5, 1:5]

#Create a data frame with PC number and percentage of variance

#Note: The percentage of variance is calculated as the squared singular value
#of each PC divided by the sum of squared singular values, multiplied by 100.

variance_table <- data.frame(
  PC = 1:length(pca$sdev),
  Variance_Percentage = pca$sdev^2 / sum(pca$sdev^2) * 100,
  cumulative_percentage = cumsum(pca$sdev^2 / sum(pca$sdev^2) * 100))

head(variance_table)

#Elbow (Scree plot)

scree <- ggplot(variance_table, aes(x = PC, y = Variance_Percentage)) +
  geom_bar(stat = 'identity', fill = "navyblue",  position = 'dodge') +
  labs(title = 'Scree plot',
       subtitle = '  ',
       x = 'Principal Components',
       y = 'Variance percentage') +
#  scale_x_discrete() +
  theme_minimal()

#Vis
scree

#Ahora hacemos la seleccion de PCs, la verdadera "reduccion de dimensionalidad"

#Table with the PCs explaining the 95% of data

variance_table_95 <- variance_table %>% 
  filter(cumulative_percentage <= 96)

#Plotting PCs explaining the 95% of data

variance_95 <- ggplot(variance_table_95,
                      aes(x = PC, y = Variance_Percentage)) +
  geom_bar(stat = 'identity', fill = 'navyblue', position = 'dodge') +
  labs(title = 'Scree plot',
       subtitle = 'for the PCs explaining the 95% of data',
       x = 'Principal Components',
       y = 'Variance percentage') +
  theme_minimal()

#Vis
variance_95

#Grafiquemos el score plot!

#Plot PC1 and PC2

PC1_PC2 <- pca_df %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point()

#Vis
PC1_PC2

#Ahora vamos a embellecerlo

#Colorear segun la metadata
PC1_PC2 <- PC1_PC2 +
  aes(colour = metadata$source_name_ch1) +
  theme_minimal()

#Vis
PC1_PC2

#Agregar titulos y subtitulos

PC1_PC2 <- PC1_PC2 +
  geom_text(mapping = aes(label = IndividualID)) +          #Agregamos texto
  labs(title = "PCA scoreplot",                           #Agregamos titulos 
       x = paste("PC1 (", paste(variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (", sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)")) #nos gusta mas usar sprintf

#Vis
PC1_PC2

#Un plus a la visualizacion, colocar elipses para ver los grupos bien definidos
PC1_PC2 <- PC1_PC2 + 
  stat_ellipse(geom = "polygon", aes(fill = metadata$source_name_ch1), alpha = 0.2, show.legend = FALSE)

#Vis
PC1_PC2

#Guardar plots

ggsave(filename = "PCA_elipse.pdf", 
       plot = PC1_PC2, 
       height = 20, 
       width = 20, 
       device = "pdf",  
       units = "cm")

#FIN
