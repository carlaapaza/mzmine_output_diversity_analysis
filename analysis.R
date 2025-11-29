###############################
# Análisis de VOCs            #
###############################

library(RColorBrewer)
library(viridis)
library(tidyverse)
library(vegan)
library(dendextend)
library(ggplot2)
library(stringr)

# Archivo de entrada (CSV) y carpeta de salida para figuras
input_file  <- "mzmine_output.csv" 
fig_dir     <- "figures"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

df <- read.csv(input_file, check.names = FALSE)
area_cols <- grep("\\.mzXML:area$", colnames(df), value = TRUE)
message("Número de columnas de área: ", length(area_cols))

int_mat <- df[, area_cols]
int_mat <- as.matrix(sapply(int_mat, as.numeric))
rownames(int_mat) <- paste0("feat_", seq_len(nrow(int_mat)))

# METADATA ----
samples <- sub("^datafile:", "", area_cols) 
samples <- sub("\\.mzXML:area$", "", samples) 
colnames(int_mat) <- samples

groups <- sub("^[^_]+_", "", samples) # el grupo está después del primer "_" en el nombre de la muestra

metadata <- data.frame(
  sample = samples,
  group  = groups,
  stringsAsFactors = FALSE
)

# TRANSFORMACIONES DE LA MATRIZ DE INTENSIDADES ----

int_mat[is.na(int_mat)] <- 0    
int_log   <- log2(int_mat + 1)  # Transformación log2 para comprimir rangos
int_log_t <- t(int_log) # filas = muestras, columnas = features

# Escalado por feature
int_scaled   <- scale(int_log) 
int_scaled_t <- t(int_scaled)

# PCA ----
pca <- prcomp(int_scaled_t, center = TRUE, scale. = FALSE)

var_exp <- summary(pca)$importance[2, 1:5]
message("Varianza explicada por los primeros PCs:")
print(var_exp)

grp        <- metadata$group
grp_levels <- unique(grp)

cols <- setNames(
  viridis(length(grp_levels), option = "viridis"),
  grp_levels
)

png(file.path(fig_dir, "pca_pc1_pc2.png"),
    width = 2000, height = 2000, res = 300)
par(mar = c(5, 5, 4, 10))
plot(pca$x[, 1], pca$x[, 2],
     col = cols[grp],
     pch = 19,
     cex = 1.5,
     xlab = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "%)"),
     ylab = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "%)"))
par(xpd = TRUE)
legend("right",
       inset  = -0.25,
       legend = grp_levels,
       col    = cols,
       pch    = 19,
       cex    = 0.8,
       title  = "Group",
       bty    = "n")
dev.off()

#######################
# DIVERSIDAD ALFA ----
#######################

# riqueza (# metabolitos por muestra)
richness <- rowSums(int_log_t > 0)

# Shannon (sensibilidad a abundancias intermedias)
alpha_shannon <- diversity(int_log_t, index = "shannon")

# Simpson (dominancia; su inverso refleja diversidad efectiva)
alpha_simpson <- diversity(int_log_t, index = "simpson")

# Evenness (uniformidad)
evenness <- alpha_shannon / log(richness)

# Tabla final diversidad alfa
alpha_df <- data.frame(
  sample  = rownames(int_log_t),
  group   = metadata$group,
  richness = richness,
  shannon  = alpha_shannon,
  simpson  = alpha_simpson,
  evenness = evenness
)

# boxplots
alpha_df_long <- alpha_df %>%
  pivot_longer(
    cols      = c(richness, shannon, simpson, evenness),
    names_to  = "metric",
    values_to = "value"
  )

alpha_df_long$metric <- str_to_title(alpha_df_long$metric)

p_alpha <- ggplot(alpha_df_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.alpha = 0.3) +
  scale_fill_viridis_d() +
  facet_wrap(~metric, scales = "free_y") +
  theme_bw(base_size = 16) +
  xlab("") +
  ylab("Valor") +
  theme(
    legend.position    = "none",
    strip.background   = element_rect(fill = "gray90", color = NA),
    strip.text         = element_text(),
    axis.text.x        = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(fig_dir, "alpha_metrics.png"),
  plot   = p_alpha,
  width  = 12,
  height = 8,
  dpi    = 300
)

#######################
# DIVERSIDAD BETA ----
#######################

# Bray-Curtis (para PCoA e intra/inter grupo)
beta_bray <- vegdist(int_log_t, method = "bray")

# Distancia basada en correlación (para PERMANOVA y betadisper)
dist_corr <- as.dist(1 - cor(int_log, method = "pearson"))

## 1. Dendrograma coloreado por grupo ----

hc   <- hclust(beta_bray, method = "ward.D2")
dend <- as.dendrogram(hc)

group_colors <- setNames(
  viridis(length(unique(metadata$group))),
  unique(metadata$group)
)

ordered_samples <- labels(dend)

label_cols <- group_colors[
  metadata$group[match(ordered_samples, metadata$sample)]
]

dend <- dend %>%
  set("labels_col", label_cols) %>%
  set("labels_cex", 0.7) %>%
  set("branches_lwd", 1)

png(file.path(fig_dir, "dendrogram_groups.png"),
    width = 2000, height = 2000, res = 300)
plot(dend, main = "Clustering jerárquico (Bray-Curtis, Ward.D2)")
legend("topright",
       legend = names(group_colors),
       col    = group_colors,
       pch    = 19,
       bty    = "n")
dev.off()

## 2. PERMANOVA + homogeneidad de dispersión ----

permanova_res <- adonis2(dist_corr ~ group, data = metadata, permutations = 999)
print(permanova_res)

bd <- betadisper(dist_corr, metadata$group)
print(anova(bd))

## 3. PCoA (Bray-Curtis) + polígonos (convex hulls) ----

pcoa_res <- cmdscale(beta_bray, eig = TRUE, k = 2)

pcoa_df <- data.frame(
  PC1   = pcoa_res$points[, 1],
  PC2   = pcoa_res$points[, 2],
  group = metadata$group
)

# función para convex hull
get_hull <- function(df) df[chull(df$PC1, df$PC2), ]

# hulls por grupo
hulls <- pcoa_df %>%
  group_by(group) %>%
  group_modify(~ get_hull(.x))

p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = group)) +
  geom_polygon(
    data        = hulls,
    aes(fill = group, color = group),
    alpha       = 0.2,
    linewidth   = 0.6,
    show.legend = FALSE
  ) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  theme_bw(base_size = 16) +
  labs(
    title = "PCoA (Bray-Curtis)",
    x     = paste0("PCoA1 (", round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1), "%)"),
    y     = paste0("PCoA2 (", round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1), "%)")
  ) +
  theme(
    legend.position = "right",
    panel.grid      = element_blank()
  )

ggsave(
  file.path(fig_dir, "pcoa_bray_poligonos.png"),
  p_pcoa,
  width  = 8,
  height = 6,
  dpi    = 300
)

## 4. Distancias intra vs inter-grupo (Bray-Curtis) ----

beta_mat <- as.matrix(beta_bray)

pairs <- expand.grid(
  sample1          = rownames(beta_mat),
  sample2          = rownames(beta_mat),
  KEEP.OUT.ATTRS   = FALSE,
  stringsAsFactors = FALSE
)

# Quedarnos solo con pares únicos (i < j)
pairs <- pairs[pairs$sample1 < pairs$sample2, ]

# Agregar distancia
pairs$dist <- mapply(
  function(x, y) beta_mat[x, y],
  pairs$sample1, pairs$sample2
)

# Agregar grupos
pairs$group1 <- metadata$group[match(pairs$sample1, metadata$sample)]
pairs$group2 <- metadata$group[match(pairs$sample2, metadata$sample)]

# Clasificar pares intra/inter
pairs$type <- ifelse(pairs$group1 == pairs$group2, "intra-group", "inter-group")

p_dist <- ggplot(pairs, aes(x = type, y = dist, fill = type)) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_fill_manual(
    values = c(
      "intra-group" = "#5ec962",
      "inter-group" = "#3b5285"
    )
  ) +
  theme_bw(base_size = 16) +
  ylab("Distancia Bray-Curtis") +
  xlab("") +
  ggtitle("Distancias intra vs. inter-grupo")

ggsave(
  file.path(fig_dir, "bray_distancias.png"),
  p_dist,
  width  = 6,
  height = 6,
  dpi    = 300
)

message("Análisis completado. Figuras guardadas en: ", fig_dir)
