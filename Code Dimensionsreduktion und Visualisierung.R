#---------------------------------------------
# 1. Pakete installieren und laden
#---------------------------------------------
install.packages(c("readr","tidyverse","Rtsne","uwot"))
library(readr)       # CSV/TSV einlesen
library(tidyverse)   # Data-Wrangling & ggplot2
library(Rtsne)       # t-SNE
library(uwot)        # UMAP
#---------------------------------------------
# 2. Datensatz einlesen
#---------------------------------------------
# Angenommen: Deine Datei heißt "data.csv" und hat die Spalten
# gene_id, sample, count, cell, dex, avgLength
pfad <- file.path("C:", "Users", "fleli", "Downloads", "data_long.csv")
df   <- read_csv(pfad, col_names = TRUE)

# Kontrollblick
glimpse(df)

#---------------------------------------------
# 3. Matrix gene × sample aufbauen
#---------------------------------------------
# Wir wollen für jede Probe ("sample") die Counts aller Gene als Features
mat_wide <- df %>%
  select(gene_id, sample, count) %>%
  pivot_wider(
    names_from  = gene_id,
    values_from = count,
    values_fill = 0
  )

# Metadaten (für Farben / Formen beim Plotten)
meta <- df %>%
  select(sample, cell, dex, avgLength) %>%
  distinct()

# Zeilenname = sample
mat_counts <- mat_wide %>%
  column_to_rownames("sample") %>%
  as.matrix()

#---------------------------------------------
# 4. Normalisierung / Transformation
#---------------------------------------------
# z. B. log1p-Transformation und Spalten‐Z-Score
#---------------------------------------------
# 4. Normalisierung / Transformation
#---------------------------------------------
mat_log <- log1p(mat_counts)                # log(1 + count)

# –––––– neu: Gene mit Varianz 0 entfernen ––––––
gene_var  <- apply(mat_log, 2, var)         # Varianz pro Spalte berechnen
mat_log   <- mat_log[, gene_var > 0]        # nur Gene mit Varianz > 0 behalten

# dann skalieren
mat_scaled <- scale(mat_log, center = TRUE, scale = TRUE)

#---------------------------------------------
# 5. Dimensionsreduktion
#---------------------------------------------

# 5.1 PCA
pca_res <- prcomp(mat_scaled, rank. = 2)
pca_coords <- as_tibble(pca_res$x[,1:2], rownames = "sample") %>%
  rename(PC1 = PC1, PC2 = PC2) %>%
  left_join(meta, by = "sample")

# 5.2 t-SNE
#---------------------------------------------
# … Dein Code bis einschließlich mat_scaled …
#---------------------------------------------
# mat_scaled hast Du bereits definiert

#---------------------------------------------
# 5.2 t-SNE mit automatischer Perplexity
#---------------------------------------------
# Anzahl Proben
n_samples   <- nrow(mat_scaled)

# maximal zulässige Perplexity (≈ (n_samples-1)/3)
max_perp    <- floor((n_samples - 1) / 3)

# wählen: Standard 30, ansonsten maximal möglichen Wert
chosen_perp <- if (max_perp >= 5) 30 else max_perp

# Ausgabe in Konsole, was wir verwenden
message("t-SNE mit perplexity = ", chosen_perp, " (n_samples = ", n_samples, ")")

set.seed(42)
tsne_res <- Rtsne(
  mat_scaled,
  dims       = 2,
  perplexity = chosen_perp,
  verbose    = TRUE
)

tsne_coords <- as_tibble(tsne_res$Y, .name_repair = "minimal") %>%
  set_names(c("TSNE1","TSNE2")) %>%
  mutate(sample = rownames(mat_scaled)) %>%
  left_join(meta, by = "sample")


# 5.3 UMAP
#---------------------------------------------
# 5.3 UMAP mit dynamischem n_neighbors
#---------------------------------------------
# Bestimme Anzahl Samples
n_samples        <- nrow(mat_scaled)

# Maximal möglicher n_neighbors (muss < n_samples)
max_n_neighbors  <- n_samples - 1

# Standardwert 15, aber nicht größer als max_n_neighbors
chosen_n_neighbors <- if (max_n_neighbors >= 15) 15 else max_n_neighbors

message("UMAP mit n_neighbors = ", chosen_n_neighbors,
        " (n_samples = ", n_samples, ")")

# UMAP ausführen
umap_res <- umap(
  mat_scaled,
  n_neighbors = chosen_n_neighbors,
  n_components = 2,
  min_dist     = 0.1,
  metric       = "cosine"
)

# Koordinaten aufbereiten
umap_coords <- as_tibble(umap_res, .name_repair = "minimal") %>%
  set_names(c("UMAP1","UMAP2")) %>%
  mutate(sample = rownames(mat_scaled)) %>%
  left_join(meta, by = "sample")


#---------------------------------------------
# 6. Visualisierung mit ggplot2
#---------------------------------------------

library(ggplot2)

# PCA-Plot
ggplot(pca_coords, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size = 3) +
  labs(title = "PCA: samples in 2D", x = "PC1", y = "PC2") +
  theme_minimal()

# t-SNE-Plot
ggplot(tsne_coords, aes(x = TSNE1, y = TSNE2, color = dex, shape = cell)) +
  geom_point(size = 3) +
  labs(title = "t-SNE: samples in 2D", x = "t-SNE1", y = "t-SNE2") +
  theme_minimal()

# UMAP-Plot
ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = dex, shape = cell)) +
  geom_point(size = 3) +
  labs(title = "UMAP: samples in 2D", x = "UMAP1", y = "UMAP2") +
  theme_minimal()
