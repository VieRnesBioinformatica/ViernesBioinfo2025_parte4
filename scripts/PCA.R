#How to PCA 
#Script to inspect and pre-process RNA-seq data counts 
#https://github.com/paulinapglz99

#I used the reference to build a PCA https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_03_Step_By_Step_PCA.pdf

#Instalación de paquetes necesarios --- ---
# Solo debe ejecutarse la primera vez si no tiene los paquetes instalados.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Instalar los paquetes GEOquery y Biobase
BiocManager::install(c("GEOquery", "Biobase"), force = TRUE)

pacman::p_load("dplyr",         #Para manejar los datos
               "ggplot2",       #Para graficar
               "stringr",       #Para manejar los nombres de genes
               "gridExtra",     #Para hacer un panel 
               "GEOquery",      #Para descargar los datos
               "Biobase",       #Para descargar los datos
               "vroom", 
               "ggfortify",
               "tibble")         #Para descargar los datos

#Descarga del conjunto de datos de expresion de GEO (GSE119834) (Pediatric Glioblastoma)

expresion <- vroom::vroom("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE119834&format=file&file=GSE119834_norm_counts_TPM_GRCh38.p13_NCBI.tsv.gz")
dim(expresion)

#Descarga el objeto GEO para obtener la metadata
gse <- getGEO("GSE119834", GSEMatrix = TRUE, AnnotGPL = FALSE)

metadata <- pData(gse[[1]])
