
# Viernes de Bioinformática -----------------------------------------------
# Sesión 7
# Instructor: LBt. Josué Guzmán Linares
#             Laboratorio Internacional EPIGEN


# Instalación y carga de paquetes necesarios ------------------------------

if (!require("pacman")) {
  install.packages("pacman")
  library(pacman)
} 


p_load("vroom",
       "dplyr",
       "circlize",
       "ComplexHeatmap",
       "tidyHeatmap")



# Heatmap con tidyHeatmap -------------------------------------------------

# tidyHeatmap está basado en el paquete ComplexHeatmap
# usa una estructura de datos basado en tidyverse. 
# Ya no es necesario convertir conjuntos de datos ordenados en matrices

# En este tutorial utilizaremos datos de espectrometría de masas
# publicados por Xia et al, Molecular neurodegeneration, 2022
# https://molecularneurodegeneration.biomedcentral.com/articles/10.1186/s13024-022-00547-7
# y reproduciremos las figuras 5a y 5b. 


# Manipulación de datos con dplyr -----------------------------------------

# Necesitaremos 2 tablas de datos

# 1. Datos de abundancias normalizadas: 
# La tabla complementaria 7 del artículo es un archivo Excel 
# que contiene las abundancias normalizadas de 225 analitos medidos en 30 muestras.
# Usaremos datos ya procesados de esa tabla complementaria 
normalized_data <- vroom("https://raw.githubusercontent.com/tomsing1/blog/refs/heads/main/posts/tidyHeatmap/abundances.csv")

# 2. Estadísticas de análisis de abundancia diferencial: 
# La tabla complementaria 8 es un archivo Excel
# con los resultados del análisis de abundancia diferencial 
# (realizado con el paquete limma de Bioconductor), 
# proporciona log2fc, valores p y valores p ajustados (FDR) para múltiples comparaciones de interés.
# Usaremos datos ya procesados de esa tabla complementaria
diff_abundance_data <- vroom("https://raw.githubusercontent.com/tomsing1/blog/refs/heads/main/posts/tidyHeatmap/statistics.csv")


# Xia et al solo incluyeron analitos con un log2fc > 0.2
# y valores p ajustados < 0,1 en la figura 5. 
# Identifiquemos los nombres de los analitos significativos en nuestro data.frame 
# de estadísticas aplicando estos umbrales. 
# A continuación, ordenamos los analitos por su log2fc y extraemos la columna component_name.
significant <- diff_abundance_data %>% 
  dplyr::filter(abs(diff_abundance_data$logFC) > 0.2, adj.P.Val < 0.1) %>% 
  dplyr::arrange(desc(logFC)) %>% 
  dplyr::pull(component_name)

head(significant)

# Conservamos 117 analitos significativos, por ejemplo, GM3(d36:1), TG(20:4_36:3) y espermina, 
# que presentan un aumento más pronunciado en APP-SAA, metoxi+, 
# en comparación con las muestras WT, metoxi-.


# A continuación, subdividimos la tabla de abundancia (normalized_data) 
# para conservar solo estos analitos significativos:
filtered_data <- normalized_data %>% 
  dplyr::filter(component_name %in% significant)


# Ahora tenemos los datos en una tabla ordenada con 3510 filas y 13 columnas.
# Para controlar el orden de las variables categóricas (por ejemplo, genotipo o lote),
# volvamoslas factores. 
# Para que coincida con el orden de los grupos de la figura 5, 
# recodificamos los niveles de la variable de grupo, 
# de modo que el grupo WT, Methoxy- aparezca en primer lugar:
filtered_data$group <-  dplyr::recode_factor(
  filtered_data$group, 
  "WT, Methoxy_Neg" = "WT, Methoxy-", 
  "APP_SAA_Hom, Methoxy_Neg" = "APP-SAA, Methoxy-",
  "APP_SAA_Hom, Methoxy_Pos" = "APP-SAA, Methoxy+"
)



# Creación de heatmaps con tidyHeatmap ------------------------------------
# El paquete tidyHeatmap proporciona funciones basadas en otro paquete llamado ComplexHeatmap
# Admite flujos de manipulación de datos basados en tidyverse, 
# además de que no necesita la conversión de tablas a matrices

# Crear un primer mapa de calor es tan fácil como pasar nuestra tabla filtered_data
# a la función tidyHeatmap::heatmap() 
# y especificar qué columnas contienen
# identificadores de fila, identificadores de columna y valores:
hm1 <- tidyHeatmap::heatmap(.data = filtered_data, 
                     .row = component_name, 
                     .column = sample_id, 
                     .value = abundance, 
                     scale = "row")

hm1

# Las abundancias normalizadas son difíciles de interpretar
# (porque los analitos se normalizaron según diferentes estándares internos
# y no se corresponden con concentraciones absolutas).
# Afortunadamente, lo que más nos interesa son los cambios relativos entre muestras 
# (y grupos de muestras). 
# Al establecer scale = «row», los valores de cada analito se convierten en z scores, 
# y la escala de colores indica la variación en torno a la media de cada fila.


# tidyHeatmap admite argumentos adicionales del paquete ComplexHeatmap::Heatmap(), 
# por ejemplo, especificar otro método para agrupar filas 
# mediante el argumento clustering_method_rows, 
# suprimir la agrupación de columnas mediante cluster_columns=FALSE, 
# o definir la paleta de colores que se va a utilizar. 

abundance_colors <- c("navy", "white", "firebrick")

hm2 <- heatmap(.data = filtered_data,
        .row = component_name, 
        .column = sample_id, 
        .value = abundance, 
        column_title = "Xia et al: targeted LC/MS data",
        row_title = "Analyte",
        scale = "row",
        cluster_columns = FALSE,
        clustering_method_rows = "ward.D", 
        palette_value = abundance_colors
)

hm2


# También podemos introducir divisiones visuales entre las columnas, 
# por ejemplo, separando los tres genotipos entre sí, 
# y asignar un color diferente a cada uno:

# Primero necesitamos usar la función group_by 
# y agrupar la columna group de filtered_data
grouped_data <- filtered_data %>% 
  group_by(group)


hm3 <- heatmap(.data = grouped_data,
          .row = component_name, 
          .column = sample_id, 
          .value = abundance, 
          column_title = "Xia et al: targeted LC/MS data",
          row_title = "Analyte",
          scale = "row",
          cluster_columns = FALSE,
          clustering_method_rows = "ward.D",
          palette_value = c("navy", "white", "firebrick"),
          palette_grouping = list(
            c("cyan4", "orange2", "darkred")
          )
  )

hm3


# Añadir anotaciones

# La función annotation_tile() nos permite añadir anotaciones adicionales 
# para identificar columnas de los lotes 
# o correspondientes al género masculino o femenino 
# en una barra de anotaciones en la parte superior del heatmap, 
# o filas con valores de diferentes paneles (lípidos, metabolitos) 
# en una franja de anotaciones a la izquierda del heatmap.

hm4 <- heatmap(.data = grouped_data, 
               .row = component_name, 
               .column = sample_id, 
               .value = abundance, 
               column_title = "Xia et al: targeted LC/MS data",
               row_title = "Analyte",
               scale = "row",
               cluster_columns = FALSE,
               clustering_method_rows = "ward.D",
               palette_value = c("navy", "white", "firebrick"),
               palette_grouping = list(
                 c("cyan4", "orange2", "darkred")
               )) %>% 
  annotation_tile(batch, palette = c("darkgreen", "gold", "darkorchid") ) %>% 
  annotation_tile(sex, palette = c("limegreen", "maroon")) %>% 
  annotation_tile(panel, palette = c("sienna", "skyblue")) %>% 
  annotation_point(cell_number)

# Algunas anotaciones son cuantitativas y pueden comunicarse mejor en un gráfico. 
# En esta nueva versión, mostramos el número de células 
# que se analizaron en cada muestra en la parte superior
# de nuestro mapa de calor con la función add_point():

hm4

# La función annotation_tile() no requiere que especifiquemos 
# si estamos añadiendo anotaciones de fila o de columna. 
# En su lugar, lo deduce automáticamente:

# ¿Cada valor utilizado para definir una fila del mapa de calor 
# (por ejemplo, component_name) se corresponde con un único valor 
# en la columna seleccionada del data.frame?
#     Si es así, cree una anotación de fila.

# ¿Cada valor utilizado para definir una columna del mapa de calor
# (por ejemplo, sample_id) se corresponde con un único valor
# en la columna seleccionada del data.frame?
#     Si es así, crea una anotación de columna.

# Si ninguna de las dos condiciones anteriores se cumple, se genera un error.



# Ordenar filas manualmente
# Para la figura 5, Xia et al ordenaron manualmente las filas, 
# colocando en la parte superior los analitos que mostraban los mayores valores log2fc
# y en la parte inferior los que presentaban las mayores disminuciones
# (es decir, cambios log2fc negativos).
# Para reproducir esto, desactivamos la agrupación de las filas 
# y cambiamos el orden del factor component_name al vector de caracteres significativo, 
# que está ordenado de mayor a menor cambio, con la función dplyr::mutate().

# reordenamos
ordered_data <- filtered_data %>% 
  dplyr::mutate(
    component_name = factor(component_name, levels = significant)
  ) %>% 
  group_by(group, panel)

# En este paso, cambiamos el comando group_by: 
# agrupamos tanto por el grupo como por las anotaciones del panel, 
# para dividir las columnas y filas de nuestro mapa de calor en paneles separados, 


# Debido a que la clusterización de filas se desactivó
# se debe especificar los colores de agrupamiento por filas 
# en la funcion palette_grouping 
hm5 <- heatmap(.data = ordered_data, 
               .row = component_name, 
               .column = sample_id, 
               .value = abundance, 
               column_title = "Xia et al: targeted LC/MS data",
               row_title = "Analyte",
               scale = "row",
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               palette_value = c("navy", "white", "firebrick"),
               palette_grouping = list(
                 c("sienna", "skyblue"),
                 c("cyan4", "orange2", "darkred")
               )) %>% 
  annotation_tile( batch, palette = c("darkgreen", "gold", "darkorchid") ) %>% 
  annotation_tile(sex, palette = c("limegreen", "maroon")) %>% 
  annotation_point(cell_number)


hm5



# Escalado manual de filas
# Hasta ahora, hemos escalado cada fila 
# convirtiendo las mediciones de cada analito en z scores, 
# y los colores reflejan las diferencias con respecto 
# a la media global de todas las muestras.

# En muchos experimentos, incluido este, 
# hay un grupo experimental explícito que proporciona una referencia natural. 
# Los efectos del genotipo APP-SAA suelen interpretarse 
# en relación con el grupo de control WT. 
# Por lo tanto, calculemos las abundancias restando el valor medio observado 
# en este grupo (por ejemplo, el WT, metoxi). 
# Con este nuevo cálculo, los colores ahora deberían indicar los cambios logarítmicos
# relativos a la mediana en el grupo de control.
new_data <- filtered_data %>% 
  dplyr::mutate(
    component_name = factor(component_name, 
                            levels = significant)
  ) %>% 
  dplyr::group_by(component_name) %>% 
  dplyr::mutate(
    ctrl_median = median(abundance[group == "WT, Methoxy-"], na.rm = TRUE),
    abundance = abundance - ctrl_median
  ) %>% 
  dplyr::ungroup() %>% 
  group_by(group, panel) 

# Aquí utilizamos funciones del paquete dplyr para

# group_by: agrupar el dataframe por component_name,
# de modo que podamos realizar cálculos por separado para cada analito,

# mutate: calcular la mediana de las abundancias observadas en el grupo WT, Methoxy- 
# mutate: restar la mediana de control de cada valor (para este analito).

# ungroup: desagrupar el dataframe, 
# para poder continuar con el conjunto de datos completo.


# agreguemos unos últimos argumentos
# para asemejar lo mejor posible a la figura 5a y 5b

# row_names_gp = grid::gpar(fontsize), aumenta el tamaño de letra de las filas

# eliminemos el nombre de columnas

# creemos un objeto que contenga la paleta de colores de abundance
# escalemos la barra de color de -2 a 2 con circlize::colorRamp2
colors <- circlize::colorRamp2(
  breaks = c(-2, 0, 2),
  colors = c("navy", "white", "firebrick")
)

# Cambiemos el título de la leyenda de los valores de abundance

hm6 <- tidyHeatmap::heatmap(.data = new_data, 
               .row = component_name, 
               .column = sample_id, 
               .value = abundance, 
               column_title = "Xia et al: targeted LC/MS data",
               row_title = "Analyte",
               row_names_gp = grid::gpar(fontsize = 8),
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               show_column_names = FALSE,
               palette_value = colors,
               heatmap_legend_param = list(title = "log2FC\n(abundance)") ,
               palette_grouping = list(
                 c("sienna", "skyblue4"),
                 c("cyan4", "orange2", "darkred")
               )) %>% 
  annotation_tile( batch, palette = c("darkgreen", "gold", "darkorchid") ) %>% 
  annotation_tile(sex, palette = c("limegreen", "maroon"))

hm6

# exportemos como png y svg desde el panel de plots de RStudio