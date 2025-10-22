# Instalar y cargar pacman (solo una vez)
if (!require("pacman")) install.packages("pacman")
library(pacman)

# Cargar (e instalar si falta) todas las librerías necesarias
p_load( dplyr, tidyr, ggplot2, purrr, vroom, stringr,
  miRNAtap, miRNAtap.db, miRBaseConverter, multiMiR, clusterProfiler,
  org.Mm.eg.db, DOSE, org.Hs.eg.db, igraph, ggraph
)

################################PREDICCIÓN DE GENES######################################################################

# Cargar el archivo CSV con los miRNAs
data <- read.csv("C:/Users/andre/Desktop/LIIGH-UNAM/GENES_MIRNAS_MIMAT.csv")

# Definir el directorio de salida
output_dir <- "C:/Users/andre/Desktop/LIIGH-UNAM/"

# Verificar la estructura del archivo
head(data)

# Obtener los IDs de miRNAs de la columna Geneid
genes_selected <- data$Gene.id

# Remover cualquier extensión "_(número)" de los IDs de miRNAs
genes_cleaned <- genes_selected %>%
  str_remove("_\\d+$")  # Elimina "_(número)" al final de los IDs

# Convertir MIMAT a nombres de miRNAs
converted_genes <- miRNA_AccessionToName(genes_cleaned, targetVersion = "v22") # Cambia la versión si es necesario

# Crear tabla de conversión
conversion_table <- data.frame(
  MIMAT_ID = genes_selected,  # Incluye los IDs originales
  Cleaned_MIMAT = genes_cleaned, # Incluye los IDs sin "_(número)"
  hsa_miR = converted_genes   # Nombres convertidos
)

# Mostrar la tabla de conversión
print(conversion_table)

# Crear la ruta completa para guardar el archivo
output_file <- file.path(output_dir, "1_MODIFICADO_miRNA_conversion_mimat_to_hsa.csv")

# Guardar el archivo en el directorio especificado
write.csv(conversion_table, output_file, row.names = FALSE)

cat("La tabla de conversión se ha guardado en:", output_file, "\n")

#########################################PREDICCIÓN DE GENES##########################################################################

# Extraer la lista de miARNs únicos de la columna mmu_miR
miRNAs_list <- unique(conversion_table$hsa_miR.TargetName)

# Filtrar los miARNs que no sean NA ni strings vacíos
miRNAs_list <- miRNAs_list[!is.na(miRNAs_list) & miRNAs_list != ""]

# Verificar los miARNs únicos cargados
cat("Lista final de miRNAs únicos para predicción:\n")
print(miRNAs_list)

# Obtener los genes blancos para todos los miARNs usando parámetros específicos
miRNAtap_results <- lapply(miRNAs_list, function(mir) {
  getPredictedTargets(
    mirna = mir, 
    sources = c("pictar", "diana", "targetscan", "miranda", "mirdb"), 
    species = "hsa", 
    min_src = 2,      # Requiere al menos 2 fuentes para considerar un target
    method = "geom",  # Método de agregación geométrica
    promote = TRUE, 
    synonyms = TRUE, 
    both_strands = FALSE
  )
})

# Nombrar los resultados con los nombres de los miARNs
names(miRNAtap_results) <- miRNAs_list

# Ver los resultados para el primer miARN
head(miRNAtap_results[[1]])

# Combinar todos los resultados en un único marco de datos
library(dplyr)
all_results <- bind_rows(lapply(names(miRNAtap_results), function(mir) {
  if (!is.null(miRNAtap_results[[mir]])) {  # Ignorar si no hay resultados
    data.frame(miRNA = mir, 
               GeneName = rownames(miRNAtap_results[[mir]]),  # Usar los nombres de los genes
               miRNAtap_results[[mir]],
               stringsAsFactors = FALSE)
  }
}), .id = NULL)

# Guardar el marco de datos como archivo CSV
write.csv(all_results, file = "2_MODIFICADO_miRNA_target_predictions_with_names.csv", row.names = FALSE)

# Crear una tabla con el número de genes blanco por cada miRNA
summary_table <- all_results %>%
  group_by(miRNA) %>%
  summarise(n_targets = n())

# Mostrar la tabla
print(summary_table)

# Seleccionar los genes blanco de un miRNA específico (por ejemplo, "mmu-miR-706")
miRNA_especifico <- "hsa-let-7a-5p"  # Cambia al miRNA que te interese
genes_de_mirna <- all_results %>% filter(miRNA == miRNA_especifico)

# Mostrar los genes blanco asociados con nombres reales
print(genes_de_mirna$GeneName)

#NUMERO DE GENES BLANCOS POR MIRNAS

# Gráfico de barras
library(ggplot2)



# Asegúrate de tener cargada la tabla summary_table
# summary_table <- all_results %>%
#   group_by(miRNA) %>%
#   summarise(n_targets = n())

# Gráfico de barras verticales con mejoras visuales y etiquetas ajustadas
ggplot(summary_table, aes(x = reorder(miRNA, -n_targets), y = n_targets, fill = n_targets)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +  # Barras con borde
  geom_text(aes(label = n_targets), vjust = -0.3, size = 3.2, color = "black") +  # Etiquetas encima
  scale_fill_gradient(low = "#A6CEE3", high = "#084594") +  # Gradiente azul
  theme_classic(base_size = 14) +
  labs(
    title = "Cantidad de genes blanco por miRNA",
    x = "miRNA",
    y = "Número de genes blanco",
    fill = "Genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 8),  # Rotar y ajustar texto
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "right",
    plot.margin = margin(1, 1, 1.5, 1.2, "cm")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))



# Cargar la tabla generada
miRNA_targets <- read.csv("2_MODIFICADO_miRNA_target_predictions_with_names.csv")


#############################################CONVERSION_DE_FORMATO_ENTREZ_SYMBOL#########################################################################################

# Cargar la tabla con las predicciones
df_genes <- read.csv("C:/Users/andre/Desktop/LIIGH-UNAM/2_MODIFICADO_miRNA_target_predictions_with_names.csv")

# Convertir los ENTREZ ID a character (necesario para mapIds)
df_genes$GeneName <- as.character(df_genes$GeneName)

# Convertir los ENTREZ ID a SYMBOL usando org.Hs.eg.db
library(org.Hs.eg.db)
library(AnnotationDbi)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = df_genes$GeneName,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Agregar los símbolos al dataframe
df_genes$GeneSymbol <- gene_symbols

# Verificar que la conversión funcionó
head(df_genes[c("GeneName", "GeneSymbol")])

# (Opcional) eliminar genes sin símbolo
df_genes <- df_genes[!is.na(df_genes$GeneSymbol), ]

# 🔹 Guardar el archivo final con la columna GeneSymbol
write.csv( df_genes, "C:/Users/andre/Desktop/LIIGH-UNAM/3_FINAL_miRNA_target_predictions_with_SYMBOL.csv", row.names = FALSE)

cat("✅ Archivo final guardado con los nombres de genes en formato SYMBOL.\n") 



#####################################HASTA AQUI#################################################33




##########Automatizado 
# Definir la ruta donde quieres guardar las imágenes
directorio_imagenes <- "C:/Users/andre/Desktop/LIIGH-UNAM/"

# Definir las categorías GO
categorias_GO <- c("BP", "CC", "MF")

# Realizar análisis de enriquecimiento GO (usando ENTREZ ID como antes)
resultados_enriquecimiento <- lapply(categorias_GO, function(ontologia) {
  enrichGO(
    gene = df_genes$GeneName,  # Usamos los ENTREZ ID para el análisis
    OrgDb = org.Hs.eg.db, 
    keyType = "ENTREZID", 
    ont = ontologia, 
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )
})

# Asignar nombres a los resultados
names(resultados_enriquecimiento) <- categorias_GO

# Generar y guardar gráficos para cada categoría (igual que antes)
for (ontologia in categorias_GO) {
  # Obtener resultados para la ontología actual
  enrich_GO <- resultados_enriquecimiento[[ontologia]]
  
  # Gráfico de barras
  bar_plot <- ggplot(enrich_GO, aes(x = reorder(Description, -FoldEnrichment), y = FoldEnrichment)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste("GO Term Enrichment -", ontologia), 
         x = "GO Terms", 
         y = "Fold Enrichment") +
    theme_minimal()
  
  # Guardar gráfico de barras
  ggsave(paste0(directorio_imagenes, "GO_Enrichment_", ontologia, "_BarPlot.png"), plot = bar_plot, width = 8, height = 6)
  
  # Dotplot
  dot_plot <- dotplot(enrich_GO, showCategory = 10) + 
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1), 
      axis.text = element_text(size = 10, color = "black"),             
      axis.title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = paste("Análisis de Enriquecimiento Funcional", ontologia),
      x = "Número de Genes",
      y = "Términos de GO"
    ) +
    scale_color_gradient(low = "blue", high = "red")
  
  # Guardar dotplot
  ggsave(paste0(directorio_imagenes, "GO_Enrichment_", ontologia, "_DotPlot.png"), plot = dot_plot, width = 8, height = 6)
}
