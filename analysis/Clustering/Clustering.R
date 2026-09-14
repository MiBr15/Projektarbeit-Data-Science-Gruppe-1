# ============================================================
# Clustering-Analyse des Airway RNA-Seq Datensatzes
# ============================================================

library(readr)
library(dplyr)
library(tibble)
library(pheatmap)
library(cluster)

# ------------------------------------------------------------
# 1. Daten einlesen
# ------------------------------------------------------------

vst_data <- read_csv(
  "data/VST_transformierte_Daten.csv",
  show_col_types = FALSE
)

metadata <- read_csv(
  "data/Probeninformationen_Normalisierung.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------
# 2. Expressionsmatrix und Metadaten vorbereiten
# ------------------------------------------------------------

vst_matrix <- vst_data %>%
  column_to_rownames("gene_id") %>%
  as.matrix()

metadata <- metadata %>%
  slice(match(colnames(vst_matrix), sample))

stopifnot(
  identical(metadata$sample, colnames(vst_matrix))
)

metadata <- metadata %>%
  mutate(
    cluster_label = paste(cell, dex, sep = "_")
  )

sample_labels <- metadata$cluster_label
names(sample_labels) <- metadata$sample

# Keine zusätzliche z-Standardisierung:
# Die Expressionswerte liegen bereits VST-transformiert auf einer gemeinsamen
# Skala vor. Eine genweise Standardisierung würde jedem Gen unabhängig von
# seiner biologischen Variabilität ein vergleichbares Gewicht geben und damit
# eine andere Fragestellung untersuchen.

# ------------------------------------------------------------
# 3. Euklidische Distanz und hierarchisches Clustering
# ------------------------------------------------------------

sample_dist <- dist(
  t(vst_matrix),
  method = "euclidean"
)

hc_complete <- hclust(
  sample_dist,
  method = "complete"
)

hc_complete_labels <- hc_complete
hc_complete_labels$labels <- sample_labels[
  hc_complete_labels$labels
]

plot(
  hc_complete_labels,
  main = "Hierarchisches Clustering der Airway-Proben",
  xlab = "Zelllinie und Behandlung",
  ylab = "Euklidische Distanz",
  sub = ""
)

# Clusterzuordnung mit dem bekannten Treatment-Label vergleichen
cluster_complete <- cutree(
  hc_complete,
  k = 2
)

table_complete <- table(
  Cluster = cluster_complete,
  Treatment = metadata$dex[
    match(names(cluster_complete), metadata$sample)
  ]
)

print(table_complete)

# ------------------------------------------------------------
# 4. Robustheitsprüfung verschiedener Linkage-Verfahren
# ------------------------------------------------------------

hc_single <- hclust(
  sample_dist,
  method = "single"
)

hc_average <- hclust(
  sample_dist,
  method = "average"
)

hc_ward <- hclust(
  sample_dist,
  method = "ward.D2"
)

# Prüft, ob k = 2 exakt nach Treatment trennt
is_treatment_separation <- function(cluster_assignment, metadata) {

  dex_ordered <- metadata$dex[
    match(names(cluster_assignment), metadata$sample)
  ]

  trt_cluster <- unique(
    cluster_assignment[dex_ordered == "trt"]
  )

  untrt_cluster <- unique(
    cluster_assignment[dex_ordered == "untrt"]
  )

  length(trt_cluster) == 1 &&
    length(untrt_cluster) == 1 &&
    trt_cluster != untrt_cluster
}

linkage_results <- tibble(
  linkage = c(
    "Single",
    "Average",
    "Complete",
    "Ward.D2"
  ),
  treatment_separation = c(
    is_treatment_separation(
      cutree(hc_single, k = 2),
      metadata
    ),
    is_treatment_separation(
      cutree(hc_average, k = 2),
      metadata
    ),
    is_treatment_separation(
      cutree(hc_complete, k = 2),
      metadata
    ),
    is_treatment_separation(
      cutree(hc_ward, k = 2),
      metadata
    )
  )
)

print(linkage_results)

# ------------------------------------------------------------
# 5. Feature-Auswahl anhand der Genvarianz
# ------------------------------------------------------------

gene_variance <- apply(
  vst_matrix,
  1,
  var
)

top500_genes <- names(
  sort(
    gene_variance,
    decreasing = TRUE
  )
)[1:500]

vst_top500 <- vst_matrix[
  top500_genes,
  ,
  drop = FALSE
]

sample_dist_top500 <- dist(
  t(vst_top500),
  method = "euclidean"
)

hc_top500 <- hclust(
  sample_dist_top500,
  method = "complete"
)

hc_top500_labels <- hc_top500
hc_top500_labels$labels <- sample_labels[
  hc_top500_labels$labels
]

plot(
  hc_top500_labels,
  main = "Clustering der 500 variabelsten Gene",
  xlab = "Zelllinie und Behandlung",
  ylab = "Euklidische Distanz",
  sub = ""
)

# ------------------------------------------------------------
# 6. Sensitivitätsanalyse der Feature-Anzahl
# ------------------------------------------------------------

n_genes <- c(
  100,
  250,
  500,
  1000,
  5000
)

feature_results <- lapply(
  n_genes,
  function(n) {

    selected_genes <- names(
      sort(
        gene_variance,
        decreasing = TRUE
      )
    )[1:n]

    matrix_subset <- vst_matrix[
      selected_genes,
      ,
      drop = FALSE
    ]

    dist_subset <- dist(
      t(matrix_subset),
      method = "euclidean"
    )

    hc_subset <- hclust(
      dist_subset,
      method = "complete"
    )

    clusters <- cutree(
      hc_subset,
      k = 2
    )

    tibble(
      genes = n,
      treatment_separation = is_treatment_separation(
        clusters,
        metadata
      )
    )
  }
) %>%
  bind_rows()

print(feature_results)

# ------------------------------------------------------------
# 7. Sample-Distance-Heatmap
# ------------------------------------------------------------

sample_dist_matrix <- as.matrix(
  sample_dist
)

rownames(sample_dist_matrix) <- sample_labels[
  rownames(sample_dist_matrix)
]

colnames(sample_dist_matrix) <- sample_labels[
  colnames(sample_dist_matrix)
]

pheatmap(
  sample_dist_matrix,
  cluster_rows = hc_complete_labels,
  cluster_cols = hc_complete_labels,
  main = "Sample-Distanzen der Airway-Proben"
)

# ------------------------------------------------------------
# 8. k-Means als unabhängiger Clustering-Ansatz
# ------------------------------------------------------------

# k = 2 wird als gezielter Robustheitscheck gewählt, da zwei
# Behandlungszustände vorliegen und das hierarchische Clustering zwei
# dominante Hauptcluster ergeben hat.

set.seed(2026)

kmeans_result <- kmeans(
  t(vst_matrix),
  centers = 2,
  nstart = 25
)

kmeans_table <- table(
  Cluster = kmeans_result$cluster,
  Treatment = metadata$dex
)

print(kmeans_table)

# ------------------------------------------------------------
# 9. Silhouettenanalyse des k-Means-Ergebnisses
# ------------------------------------------------------------

silhouette_result <- silhouette(
  kmeans_result$cluster,
  sample_dist
)

plot(
  silhouette_result,
  main = "Silhouette-Plot des k-Means-Clusterings"
)

mean_silhouette <- mean(
  silhouette_result[, "sil_width"]
)

print(mean_silhouette)
