# 0. (falls nötig) Pakete installieren
# install.packages(c("readr","dplyr","ggplot2","Rtsne","uwot"))

# 1. Pakete laden
library(readr)     # CSV einlesen
library(dplyr)     # Data-Wrangling
library(ggplot2)   # Plotting
library(Rtsne)     # t-SNE
library(uwot)      # UMAP

# 2. VST‐Daten importieren
pfad <- file.path("VST_transformierte_Daten.csv")
vst_df <- read_csv(pfad, col_names = TRUE)
glimpse(vst_df)    # prüfen, dass gene_id + SRR-Spalten da sind

# 3. Matrix bauen
#    – Zeilen: samples (SRR…), Spalten: GeneIDs
mat <- vst_df %>%
  column_to_rownames("gene_id") %>%  # gene_id → rownames
  as.matrix()                        # Datenmatrix (Gene × Sample)
mat_t <- t(mat)                      # transpose → Sample × Gene

# 4. Features skalieren
#    – jede Gen-Spalte auf Mittel 0 und SD 1 bringen
mat_scaled <- scale(mat_t)

# 5. Dimensionsreduktion

## 5.1 PCA
pca_res <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)
pca_df  <- as_tibble(pca_res$x[,1:2], rownames = "sample") %>%
  set_names("sample","PC1","PC2")

## 5.2 t-SNE (Perplexity automatisch anpassen)
n_samp   <- nrow(mat_scaled)
max_p     <- floor((n_samp - 1) / 3)
perplex  <- if (max_p >= 5) 30 else max_p
message("Verwende t-SNE perplexity = ", perplex)
set.seed(42)
tsne_res <- Rtsne(mat_scaled, dims = 2, perplexity = perplex, verbose = TRUE)
tsne_df  <- as_tibble(tsne_res$Y, .name_repair = "minimal") %>%
  set_names("TSNE1","TSNE2") %>%
  mutate(sample = rownames(mat_scaled))

## 5.3 UMAP (n_neighbors automatisch anpassen)
n_nbrs   <- min(15, n_samp - 1)
message("Verwende UMAP n_neighbors = ", n_nbrs)
umap_res <- umap(mat_scaled, n_neighbors = n_nbrs, min_dist = 0.1)
umap_df  <- as_tibble(umap_res, .name_repair = "minimal") %>%
  set_names("UMAP1","UMAP2") %>%
  mutate(sample = rownames(mat_scaled))

# 6. Plotten

# 6.1 PCA-Plot
ggplot(pca_df, aes(PC1, PC2, label = sample)) +
  geom_point() +
  geom_text(vjust = -0.5, size = 3) +
  labs(title = "PCA-Projektion", x = "PC1", y = "PC2") +
  theme_minimal()

# 6.2 t-SNE-Plot
ggplot(tsne_df, aes(TSNE1, TSNE2, label = sample)) +
  geom_point() +
  geom_text(vjust = -0.5, size = 3) +
  labs(title = "t-SNE-Projektion", x = "t-SNE 1", y = "t-SNE 2") +
  theme_minimal()

# 6.3 UMAP-Plot
ggplot(umap_df, aes(UMAP1, UMAP2, label = sample)) +
  geom_point() +
  geom_text(vjust = -0.5, size = 3) +
  labs(title = "UMAP-Projektion", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal()
