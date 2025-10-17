# Viernes de Bioinformática -----------------------------------------------
# Sesión 7
# Instructor: LBt. Josué Guzmán Linares
#             Laboratorio Internacional EPIGEN


# Instalación y carga de paquetes necesarios ------------------------------

# Verifica si el paquete "pacman" está instalado; si no lo está, lo instala
if (!require("pacman")) {
  install.packages("pacman")  # Instala el paquete "pacman"
  library(pacman)             # Carga la librería una vez instalada
}

# Carga múltiple de paquetes con pacman::p_load()
# vroom: para lectura rápida de archivos grandes
# dplyr: para manipulación de datos
# ggplot2: para visualización de datos
p_load("vroom",
       "dplyr",
       "ggplot2")


# Descarga y lee un archivo TSV (tab-separated values) directamente desde una URL
# Contiene datos sobre muertes por drogas en Chicago
chicago_deaths_data <- vroom("https://riffomonas.org/code_club/assets/data/chicago_drug_deaths.tsv")

# Crea un heatmap básico con ggplot2
# year: eje x, age: eje y, n: intensidad de color (muertes)
hm1 <- ggplot(data = chicago_deaths_data,
              aes(x = year, y = age, fill = n)) +
  geom_tile(colour = "white") +   # Dibuja las celdas del heatmap con borde blanco
  geom_abline(intercept = c(-1951, -1970), slope = 1,  # Agrega líneas diagonales (abline)
              linetype = "longdash", linewidth = 0.25)  # Define el tipo y grosor de línea

hm1  # Muestra el gráfico hm1


# Añade una línea vertical y un texto descriptivo al heatmap
hm2 <- hm1 +
  geom_vline(xintercept = 1988.45, color = "white") +   # Línea vertical blanca en x=1988.45
  annotate(geom = "text",                               # Agrega texto anotativo al gráfico
           hjust = 0,                                   # Alineación horizontal del texto
           x = 2005.5,                                  # Posición en eje x
           y = 31,                                      # Posición en eje y
           label = "Men born from 1951 to 1970",        # Texto que se mostrará
           size = 3)                                    # Tamaño del texto

hm2  # Muestra el gráfico hm2


# Personaliza la escala de colores, ejes y etiquetas
hm3 <- hm2 +
  scale_fill_gradientn(                                 # Escala de color personalizada
    colours = c("#F4F3E8", "#FBCFB0", "#F7AC88", "#EB8772", "#B6525B",
                "#C34280", "#9A248E", "#6E1193", "#29088C", "#0D0887",
                "#0D0887", "#0C0887", "#0D0987"),
    values = c(0, 50, 100, 150, 200,
               250, 300, 350, 400, 450,
               500, 600, 735) / 735,                    # Normaliza los valores de color (0 a 1)
    labels = c(200, 400, 600),                          # Etiquetas para la leyenda
    breaks = c(200, 400, 600)                           # Puntos de ruptura en la leyenda
  ) +
  scale_x_continuous(breaks = seq(1990, 2020, 5)) +     # Eje X con saltos cada 5 años
  scale_y_continuous(breaks = seq(10, 80, 10)) +        # Eje Y con saltos de 10 años de edad
  coord_cartesian(expand = FALSE, clip = "off",          # Ajusta los límites del gráfico
                  xlim = c(1988.4, NA)) +                # Fija límite inferior en x
  labs(
    title = "Drug deaths in Chicago among Black men",    # Título del gráfico
    caption = "Source: Times/Banner analysis of N.C.H.S. mortality data", # Fuente de los datos
    fill = "Deaths per 100k",                            # Título de la leyenda
    x = "year",                                          # Etiqueta eje X
    y = "age")                                           # Etiqueta eje Y

hm3  # Muestra el gráfico hm3


# Personaliza la apariencia del tema (theme)
hm4 <- hm3 +
  theme(
    legend.position = "bottom",                          # Coloca la leyenda abajo
    legend.justification.bottom = "right",               # Justificación de la leyenda
    legend.key.height = unit(5, "pt"),                   # Altura de las cajas de leyenda
    legend.key.width = unit(24, "pt"),                   # Ancho de las cajas de leyenda
    panel.grid = element_blank(),                        # Elimina líneas de cuadrícula
    panel.background = element_blank()                   # Fondo blanco sin panel
  )

hm4  # Muestra el gráfico final


# Guarda el gráfico final como imagen PNG con resolución 300 dpi
ggsave(plot = hm4, filename = "chicago_drug_deaths.png", width = 5, height = 6, dpi = 300)