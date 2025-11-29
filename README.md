# Análisis de VOCs (GC-MS)

Este repositorio contiene el análisis de perfiles de metabolitos volátiles (VOCs) detectados por GC-MS. El script realiza:

- **Carga y limpieza de datos** desde un archivo CSV con columnas `*.mzXML:area`.
- **Transformaciones** log2 y escalado.
- **PCA** para variación global entre muestras.
- **Diversidad alfa**: riqueza, Shannon, Simpson, evenness.
- **Diversidad beta**: Bray-Curtis, correlación 1-Pearson.
- **PCoA**, **clustering jerárquico**, **convex hulls**.
- **PERMANOVA** y **betadisper**.
- **Comparación de distancias intra vs inter-grupo**.

Todas las figuras generadas se guardan automáticamente en un directorio figures/

## Requisitos

Paquetes de R:
tidyverse
vegan
viridis
dendextend
ggplot2
RColorBrewer
