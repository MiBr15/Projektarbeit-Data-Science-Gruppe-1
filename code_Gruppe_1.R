## ============================================================
# Pakete
# ============================================================

library(airway)
library(tidyverse)
library(pheatmap)
library(edgeR)
library(DESeq2)
library(apeglm)
library(Rtsne)
library(uwot)
library(cluster)

data(airway)

# ============================================================
# Data Wrangling
# ============================================================

#1. Datenstruktur verstehen:

airway

assay(airway)

colData(airway)

head(assay(airway))

dim(airway)

#2. Datenextraktion

Messwerte_roh <- assay(airway)

Probeninformationen <- as.data.frame(colData(airway))

Geninformationen <- as.data.frame(rowData(airway))

#3. Untersuchung Messwerte_roh,Probeninformationen, Geninformationen

class(Messwerte_roh)
dim(Messwerte_roh)
head(Messwerte_roh)
colnames(Messwerte_roh)
rownames(Messwerte_roh)

class(Probeninformationen)
dim(Probeninformationen)
head(Probeninformationen)
colnames(Probeninformationen)
rownames(Probeninformationen)

class(Geninformationen)
dim(Geninformationen)
head(Geninformationen)
colnames(Geninformationen)
rownames(Geninformationen)

#4. Überprüfung gleicher Sample-IDs

all(colnames(Messwerte_roh) == Probeninformationen$Run)

#5. Erstellen einer Probeninformationstabelle

Tab_Probeninformationen <- Probeninformationen[, c(
  "Run",
  "cell",
  "dex",
  "avgLength"
)]

#6. Überprüfung erstellter Tabelle

nrow(Tab_Probeninformationen)

length(unique(Tab_Probeninformationen$Run))

#7. Messwerte_roh ins Long Format bringen

Messwerte_roh_long <- as.data.frame(Messwerte_roh)

Messwerte_roh_long$gene_id <- rownames(Messwerte_roh_long)

rownames(Messwerte_roh_long) <- NULL

Messwerte_roh_long <- pivot_longer(
  Messwerte_roh_long,
  col = -gene_id,
  names_to ="sample",
  values_to = "counts"
  
)

#8. Überprüfung der Tabelle

dim(Messwerte_roh_long)

#9. Messwerte_roh_long mit Tab_Probeninformationen verbinden

Daten_roh_long <- Messwerte_roh_long %>%
  left_join(
    Tab_Probeninformationen,
    by = c("sample" = "Run")
  )

#10. Neue Tabelle überprüfen

head(Daten_roh_long)

dim(Daten_roh_long)

sum(is.na(Daten_roh_long$cell))


#11. Join von Daten_roh_long und Geninformationen

Daten_roh_long <- Daten_roh_long %>%
  left_join(
    Geninformationen,
    by = "gene_id"
  )

#12. Überprüfung des Daten_roh_long-Datensatzes

##Aufbau der Rohdatentabelle

glimpse(Daten_roh_long)

##Anzahl der Proben

n_distinct(Daten_roh_long$sample)

##Anzahl der Gene

n_distinct(Daten_roh_long$gene_id)

##Anzahl der Gene pro Probe

Daten_roh_long %>%
  count(sample)

##Verteilung der Genprobendaten auf treated und untreated

Daten_roh_long %>%
  count(dex)

##Anzahl Genprobendaten pro Zelllinie

Daten_roh_long %>%
  count(cell, dex)

##Suche nach NA-Werten

colSums(is.na(Daten_roh_long))

#13. Rohdatentabelle als csv speichern

write.csv(
  Daten_roh_long,
  "Daten_roh_long.csv",
  row.names = FALSE
)

# ============================================================
# QA
# ============================================================

# 1. Welche Daten verwendet werden

data <- Daten_roh_long


# 2. KONTROLLE DER DATENSTRUKTUR

str(data)

head(data)

dim(data)

colnames(data)


# Anzahl NA-Werte pro Variable

colSums(is.na(data))

# Gesamtzahl der NA-Werte

sum(is.na(data))


# Duplikate?

sum(
  duplicated(
    data[, c("gene_id", "sample")]
  )
)


# sample/gen anzahl
length(unique(data$sample))

length(unique(data$gene_id))

unique(data$sample)

unique(data$cell)

unique(data$dex)


# Hier startet die QA


# 3. Count-Matrix

count_data <- data %>%
  select(
    gene_id,
    sample,
    counts
  ) %>%
  pivot_wider(
    names_from = sample,
    values_from = counts
  )

count_data <- as.data.frame(count_data)

rownames(count_data) <- count_data$gene_id

count_data$gene_id <- NULL

count_matrix <- as.matrix(count_data)


# Kontrolle

dim(count_matrix)

head(count_matrix)


# 4. Library size

library_size <- colSums(count_matrix)

library_size_df <- data.frame(
  sample = names(library_size),
  total_counts = as.numeric(library_size)
)

ggplot(
  library_size_df,
  aes(
    x = sample,
    y = total_counts
  )
) +
  geom_col(
    fill = "#08306B"
  ) +
  labs(
    title = "Library Size pro Sample",
    x = "Sample",
    y = "Total Counts"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# 5. Anzahl der detektierten Gene

detected_genes <- colSums(
  count_matrix > 0
)

zero_counts <- colSums(
  count_matrix == 0
)

sample_qc <- data.frame(
  sample = colnames(count_matrix),
  library_size = colSums(count_matrix),
  detected_genes = detected_genes,
  zero_counts = zero_counts
)

sample_qc

ggplot(
  sample_qc,
  aes(
    x = sample,
    y = detected_genes
  )
) +
  geom_col(
    fill = "#08306B"
  ) +
  labs(
    title = "Anzahl detektierter Gene pro Sample",
    x = "Sample",
    y = "Anzahl Gene mit Count > 0"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


# 6. Mean-Varianz-Analyse

gene_mean <- rowMeans(count_matrix)

gene_variance <- apply(
  count_matrix,
  1,
  var
)

mean_variance <- data.frame(
  gene = rownames(count_matrix),
  mean = gene_mean,
  variance = gene_variance
)

ggplot(
  mean_variance,
  aes(
    x = mean,
    y = variance
  )
) +
  geom_point(
    alpha = 0.2,
    color = "#08306B"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Mean-Varianz-Plot",
    x = "Mean Count (log10)",
    y = "Varianz (log10)"
  ) +
  theme_minimal(base_size = 13)


# 7. Mean-SD-Plot

gene_sd <- apply(
  count_matrix,
  1,
  sd
)

mean_sd_raw <- data.frame(
  gene = rownames(count_matrix),
  mean = gene_mean,
  sd = gene_sd
)

ggplot(
  mean_sd_raw,
  aes(
    x = mean,
    y = sd
  )
) +
  geom_point(
    alpha = 0.2,
    color = "#08306B"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Mean-SD-Plot",
    x = "Mean Count (log10)",
    y = "Standardabweichung (log10)"
  ) +
  theme_minimal(base_size = 13)


# 8. Log2 Transformation und Mean-SD-Plot mit log2 Daten

log_counts <- log2(
  count_matrix + 1
)

# mean-sd mit log2

log_mean <- rowMeans(log_counts)

log_sd <- apply(
  log_counts,
  1,
  sd
)

mean_sd_log <- data.frame(
  gene = rownames(log_counts),
  mean = log_mean,
  sd = log_sd
)

ggplot(
  mean_sd_log,
  aes(
    x = mean,
    y = sd
  )
) +
  geom_point(
    alpha = 0.2,
    color = "#08306B"
  ) +
  labs(
    title = "Mean-SD-Plot nach log2-Transformation",
    x = "Mean (log2)",
    y = "Standardabweichung"
  ) +
  theme_minimal(base_size = 13)


# ============================================================
# QC
# ============================================================

# 1. AUSGABEVERZEICHNIS

output_dir <- "data"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# 2. SPALTEN CHECKEN

required_cols <- c(
  "gene_id",
  "sample",
  "counts"
)

missing_cols <- setdiff(
  required_cols,
  colnames(data)
)

if (length(missing_cols) > 0) {
  stop(
    "Fehlende Pflichtspalten: ",
    paste(missing_cols, collapse = ", ")
  )
}


# 3. DATENTYPEN KONTROLLIEREN

data <- data %>%
  mutate(
    gene_id = as.character(gene_id),
    sample = as.character(sample),
    counts = as.numeric(counts)
  )

if (any(is.na(data$counts)) || any(data$counts < 0)) {
  stop("'counts' enthält NA oder negative Werte.")
}

data <- data %>%
  mutate(
    counts = as.integer(round(counts))
  )


# 4. COUNT-MATRIX AUS DER CSV ERSTELLEN

count_data <- data %>%
  select(
    gene_id,
    sample,
    counts
  ) %>%
  pivot_wider(
    names_from = sample,
    values_from = counts
  )

count_data <- as.data.frame(count_data)

rownames(count_data) <- count_data$gene_id

count_data$gene_id <- NULL

count_matrix <- as.matrix(count_data)


# Kontrolle

cat("\nDimension der Count-Matrix: ")
print(dim(count_matrix))


# 5. edgeR DGEList ERSTELLEN

dge <- DGEList(
  counts = count_matrix
)


# 6. LOW-COUNT-FILTERUNG MIT filterByExpr()

keep <- filterByExpr(dge)

cat("\nGene vor Filterung: ", nrow(dge), "\n")
cat("Gene nach Filterung: ", sum(keep), "\n")
cat("Entfernte Gene: ", sum(!keep), "\n")

if (sum(keep) == 0) {
  stop("Kein Gen erfüllt die filterByExpr()-Kriterien.")
}

filtered_count_matrix <- count_matrix[
  keep,
  ,
  drop = FALSE
]


# 7. GEFILTERTE GENE WIEDER IN DIE LONG-TABELLE ZURÜCKFÜHREN

filtered_gene_ids <- rownames(
  filtered_count_matrix
)

data_filtered <- data %>%
  filter(
    gene_id %in% filtered_gene_ids
  )


# 8. GEFILTERTE DATEN ALS CSV SPEICHERN

output_file <- file.path(
  output_dir,
  "Daten_roh_long_filtered.csv"
)

write_csv(
  data_filtered,
  output_file,
  na = ""
)


# 9. AUSGABE FÜR ÜBERSICHT

cat("\nLow-Count-Filterung abgeschlossen.\n")

cat(
  "Gene vor Filterung: ",
  nrow(count_matrix),
  "\n"
)

cat(
  "Gene nach Filterung: ",
  nrow(filtered_count_matrix),
  "\n"
)

cat(
  "Entfernte Gene: ",
  sum(!keep),
  "\n"
)

cat(
  "Zeilen vor Filterung: ",
  nrow(data),
  "\n"
)

cat(
  "Zeilen nach Filterung: ",
  nrow(data_filtered),
  "\n"
)

cat(
  "Ausgabedatei: ",
  output_file,
  "\n"
)


# ============================================================
# Data Processing
# ============================================================


# 1. Einlesen der gefilterten Rohdaten

fil_Daten <- data_filtered

# 2. Prüfung der überreichten Daten_roh_long_filtered.csv

##Prüfung ob alle nötigen Spalten noch vorliegen

required_cols <- c(
  "gene_id",
  "sample",
  "counts",
  "cell",
  "dex"
)

missing_cols <- setdiff(required_cols, colnames(fil_Daten))

if (length(missing_cols) > 0) {
  stop(
    "Fehlende nötige Spalten: ",
    paste(missing_cols, collapse = ", ")
  )
}

## Prüfung ob NA- oder negative Werte vorliegen und ob alle Werte Integer sind

if (any(is.na(fil_Daten$counts))) {
  stop("'counts' -> NA-Werte.")
}

if (any(fil_Daten$counts < 0)) {
  stop("'counts' -> negative Werte.")
}

if (any(fil_Daten$counts != round(fil_Daten$counts))) {
  stop("'counts' -> keine reinen Integer.")
}

## Prüfung der eindeutigen Definition aller Samples

Probeninformationen <- fil_Daten %>%
  distinct(sample, cell, dex)

if (nrow(Probeninformationen) != n_distinct(fil_Daten$sample)) {
  stop(
    "Min. ein Sample besitzt mehrere unterschiedliche ",
    "cell- oder dex-Zuordnungen."
  )
}

# 3. Extrahierung einer Count-Matrix aus fil_Daten

Daten_matrix <- fil_Daten %>%
  group_by(gene_id, sample) %>%
  summarise(
    counts = sum(counts),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = sample,
    values_from = counts,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(Daten_matrix) <- Daten_matrix$gene_id
Daten_matrix$gene_id <- NULL

Daten_matrix <- as.matrix(Daten_matrix)

storage.mode(Daten_matrix) <- "integer"

## Überprüfung der Matrix

dim(Daten_matrix)

head(Daten_matrix[, 1:4])

colnames(Daten_matrix)

# 4. Probeninformationen (Metadaten) zuordnen

Probeninformationen <- fil_Daten %>%
  distinct(sample, cell, dex) %>%
  as.data.frame()

rownames(Probeninformationen) <- Probeninformationen$sample

Probeninformationen <- Probeninformationen[
  colnames(Daten_matrix),
  ,
  drop = FALSE
]

## Prüfung der korrekten Daten-Metadatenzuordnung

stopifnot(
  identical(
    rownames(Probeninformationen),
    colnames(Daten_matrix)
  )
)

# 5. Beginn des Preprocessing/Normalisierens - Erstellen eines DESeq2-Objekts

Probeninformationen$cell <- factor(Probeninformationen$cell)
Probeninformationen$dex <- factor(Probeninformationen$dex)

dds <- DESeqDataSetFromMatrix(
  countData = Daten_matrix,
  colData = Probeninformationen,
  design = ~ cell + dex
)

# 6. Size-Factor-Normalisierung

dds <- estimateSizeFactors(dds)

sizeFactors(dds)

# 7. Normalisierung der Messdaten

normalisierte_counts <- counts(
  dds,
  normalized = TRUE
)

# 8. Kontrolle der normalisierten Counts

dim(normalisierte_counts)

head(normalisierte_counts[, 1:4])

data.frame(
  sample = colnames(Daten_matrix),
  size_factor = sizeFactors(dds)
)

# 9. Variance Stabilizing Transformation (VST)

vst_daten <- vst(
  dds,
  blind = TRUE
)

vst_matrix <- assay(vst_daten)

# 10. QA der Size Factors

size_factor_info <- data.frame(
  sample = colnames(Daten_matrix),
  cell = Probeninformationen$cell,
  dex = Probeninformationen$dex,
  size_factor = sizeFactors(dds)
)

print(size_factor_info)

# 11. QA der Library Size

library_größe <- colSums(Daten_matrix)

library_info <- data.frame(
  sample = names(library_größe),
  library_größe = as.numeric(library_größe),
  size_factor = sizeFactors(dds)
)

print(library_info)

# 12. QA der VST-Matrix

dim(vst_matrix)

head(vst_matrix[, 1:4])

# 13. QA des Mean-Variance-Plots

vst_mean <- rowMeans(vst_matrix)

vst_variance <- apply(
  vst_matrix,
  1,
  var
)

plot(
  vst_mean,
  vst_variance,
  log = "xy",
  pch = 16,
  cex = 0.4,
  xlab = "Mean VST expression",
  ylab = "Variance",
  main = "Mean-Variance relationship after VST"
)

# 14. Ergebnisse speichern

output_dir_data <- "data"
output_dir_images <- "images"

vst_ergebnisse <- as.data.frame(vst_matrix) %>%
  rownames_to_column("gene_id")

normalisierte_ergebnisse <- as.data.frame(normalisierte_counts) %>%
  rownames_to_column("gene_id")

normalisierungs_info <- Probeninformationen %>%
  mutate(
    size_factor = sizeFactors(dds)
  )

output_dir_data <- "data"
output_dir_images <- "images"

if (!dir.exists(output_dir_data)) {
  dir.create(output_dir_data, recursive = TRUE)
}

if (!dir.exists(output_dir_images)) {
  dir.create(output_dir_images, recursive = TRUE)
}



png(
  filename = file.path(
    output_dir_images, 
    "VST_Mittelwert-Varianz-Beziehung.png"
  ),
  width = 1600,
  height = 1200,
  res = 200
)

plot(
  vst_mean,
  vst_variance,
  log = "xy",
  pch = 16,
  cex = 0.4,
  xlab = "Mittlere VST-Expression",
  ylab = "Varianz",
  main = "Mittelwert-Varianz-Beziehung nach erfolgter VST"
)

dev.off()

write_csv(
  normalisierte_ergebnisse,
  file.path(
    output_dir_data,
    paste0("DESeq2_normalisierte_counts.csv")
  ),
  na = ""
)

write_csv(
  vst_ergebnisse,
  file.path(
    output_dir_data,
    paste0("VST_transformierte_Daten.csv")
  ),
  na = ""
)

write_csv(
  normalisierungs_info,
  file.path(
    output_dir_data,
    paste0("Probeninformationen_Normalisierung.csv")
  ),
  na = ""
)


# ============================================================
# Differentialexpressionsanalyse
# ============================================================


# 1. Daten aus csv laden
Daten_roh_long_filtered <- data_filtered

Tab_Probeninformationen <- as.data.frame(
  Tab_Probeninformationen
)

# 2. Daten_roh_long_filtered für DESeq2 vorbereiten. Umwandlung in eine wide-format count Matrix. Die Spaltennamen
# der count matrix müssen mit den Zeilennamen der Metadaten übereinstimmen
count_mat <- Daten_roh_long_filtered %>%
  select (gene_id, sample, counts) %>%
  pivot_wider (
    names_from = sample,
    values_from = counts
  ) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()


# 3. cell und dex als Faktoren definieren. Somit kann für die spätere Auswertung mit DESeq2 in Kategorien
# wie treated und untreated und in die vier unterschiedlichen Zelllinien/Spender eingeteilt werden.
Tab_Probeninformationen$cell <-
  factor(Tab_Probeninformationen$cell)

Tab_Probeninformationen$dex <-
  factor(Tab_Probeninformationen$dex)

rownames(Tab_Probeninformationen) <- Tab_Probeninformationen$Run

# Faktoren können überprüft werden
levels(Tab_Probeninformationen$cell)


# 4. Untreated als Referenzgruppe festlegen. So kann später intuitiver bestimmt werden, wie sich die Behandlung
# im Gegensatz zu keiner Behandlung auf die Genexpression auswirkt.
Tab_Probeninformationen$dex <- relevel(
  Tab_Probeninformationen$dex,
  ref = "untrt"
)

#Überprüfen ob untrt wirklich als Referenzlevel festgelegt wurde. Reihenfolge muss "untrt" "trt" sein
levels(Tab_Probeninformationen$dex)

# 5. Für DESeq2 müssen Reihenfolge der Count Matrix mit den Metadaten identisch sein. Dies wird überprüft. Falls es nicht stimmt
# bricht das Skript ab
stopifnot(
  identical(
    colnames(count_mat),
    rownames(Tab_Probeninformationen)
  )
)

# 6. DESeq2 durchführen
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = Tab_Probeninformationen,
  design = ~ cell+dex
  
)

dds <- DESeq(dds)
#Überprüfen, ob resultsNames(dds) dex_trt_vs_untrt entspricht, da dieser Koeffizient später noch verwendet wird
stopifnot(
  "dex_trt_vs_untrt" %in% resultsNames(dds)
)

# 7. Differentiale Genexpression zwischen behandelten (trt) und unbehandelten
# Proben (untrt) testen. Alpha = 0.05 beschreibt, dass Gene mit einem adjustierten p-value <0,05 signifikant sind.
res <- results(
  dds,
  contrast = c("dex","trt","untrt"),
  alpha = 0.05
)
summary(res)


# 8. Log2-FC shrinken, um unsichere oder möglicherweise überschätzte Effekte, besonders bei Genen mit wenigen Counts oder hoher Streuung, zu stabilisieren.
res_shrunk <- lfcShrink(
  dds,
  coef = "dex_trt_vs_untrt",
  type = "apeglm"
)

# 9. Ergebnisse sollen übersichtlicher gestaltet werden, weshalb Geninformationen dazugeholt werden. Zunächst wird eine Annotationstabelle erstellt.
annotation <- Daten_roh_long_filtered %>%
  select(
    gene_id,
    gene_name,
    symbol,
    gene_biotype
  ) %>%
  distinct()

# 10. DESeq2 Ergebnisse werden in Data Frame umgewandelt.
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id")

# 11. Annotation anhängen
res_df <- res_df %>%
  left_join(
    annotation,
    by = "gene_id"
  )

# 12. Nach Signifikanz sortieren. Kleinsten padj kommen nach oben.
res_df <- res_df %>%
  arrange(padj)

# 13. Eigene Tabelle für signifikante Gene (padj < 0.05) erstellen.
sig <- res_df %>%
  filter(
    !is.na(padj),
    padj < 0.05
  )

# 14. Übersichtlichere Tabelle mit signifikanten Werten und wichtigen Informationen.
sig <- sig %>%
  select(
    gene_id,
    gene_name,
    symbol,
    gene_biotype,
    baseMean,
    log2FoldChange,
    lfcSE,
    pvalue,
    padj
  )

# 15. Betrachtung nur hochregulierter und dann nur herunterregulierter Gene nach der Behandlung.
sig_up <- sig %>%
  filter(log2FoldChange > 0)

sig_down <- sig %>%
  filter(log2FoldChange < 0)

# 16. Sortieren nach dem stärksten Effekt. Welche Gene besitzen einen hohen LFC. Welche wurden besonders stark hochreguliert?
sig_up %>%
  arrange(desc(log2FoldChange)) %>%
  select(
    symbol,
    gene_name,
    baseMean,
    log2FoldChange,
    padj
  ) %>%
  head(20)

# 17. Welche Gene wurden am stärksten runterreguliert?
sig_down %>%
  arrange(log2FoldChange) %>%
  select(
    symbol,
    gene_name,
    baseMean,
    log2FoldChange,
    padj
  ) %>%
  head(20)

# 18.Ergebnistabelle als csv speichern
# Komplette annotierte Tabelle
write_csv(
  res_df,
  "data/DESeq2_Ergebnisse_alle_Gene.csv"
)
# Nur signifikante Gene
write_csv(
  sig,
  "data/DESeq2_signifikante_Gene.csv"
)

# 19. MA-Plot aus geshrinkten LFCs
plotMA(res_shrunk, ylim = c(-2, 2))


# ============================================================
# Dimensionsreduktion und Visualisierung
# ============================================================


# 1. Zeilen: Samples, Spalten: Gene
mat_t <- t(vst_matrix)

# 2. Features skalieren
mat_scaled <- scale(mat_t)

# 3. Dimensionsreduktion

plot_info <- Probeninformationen %>%
  as.data.frame() %>%
  select(
    sample,
    cell,
    dex
  )

## 3.1 PCA
pca_res <- prcomp(mat_scaled, center = FALSE, scale. = FALSE)
pca_variance <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
pca_df  <- as_tibble(pca_res$x[,1:2], rownames = "sample") %>%
  set_names("sample","PC1","PC2")
pca_df <- pca_df %>%
  left_join(
    plot_info,
    by = "sample"
  )

## 3.2 t-SNE (Perplexity automatisch anpassen)
n_samp   <- nrow(mat_scaled)
max_p     <- floor((n_samp - 1) / 3)
perplex  <- if (max_p >= 5) 30 else max_p
message("Verwende t-SNE perplexity = ", perplex)
set.seed(42)
tsne_res <- Rtsne(mat_scaled, dims = 2, perplexity = perplex, verbose = TRUE)
tsne_df  <- as_tibble(tsne_res$Y, .name_repair = "minimal") %>%
  set_names("TSNE1","TSNE2") %>%
  mutate(sample = rownames(mat_scaled))
tsne_df <- tsne_df %>%
  left_join(
    plot_info,
    by = "sample"
  )

## 3.3 UMAP (n_neighbors automatisch anpassen)
n_nbrs   <- min(15, n_samp - 1)
message("Verwende UMAP n_neighbors = ", n_nbrs)
set.seed(42)
umap_res <- umap(mat_scaled, n_neighbors = n_nbrs, min_dist = 0.1)
umap_df  <- as_tibble(umap_res, .name_repair = "minimal") %>%
  set_names("UMAP1","UMAP2") %>%
  mutate(sample = rownames(mat_scaled))
umap_df <- umap_df %>%
  left_join(
    plot_info,
    by = "sample"
  )

# 4. Plotten

# 4.1 PCA-Plot
ggplot(
  pca_df,
  aes(
    PC1,
    PC2,
    label = sample,
    color = dex,
    shape = cell
  )
) +
  geom_point(size = 3) +
  geom_text(
    vjust = -0.7,
    size = 3
  ) +
  labs(
    title = "PCA-Projektion",
    x = "PC1",
    y = "PC2",
    color = "Behandlung",
    shape = "Zelllinie"
  ) +
  theme_minimal()

# 4.2 t-SNE-Plot
ggplot(
  tsne_df,
  aes(
    TSNE1,
    TSNE2,
    label = sample,
    color = dex,
    shape = cell
  )
) +
  geom_point(size = 3) +
  geom_text(
    vjust = -0.7,
    size = 3
  ) +
  labs(
    title = "t-SNE-Projektion",
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Behandlung",
    shape = "Zelllinie"
  ) +
  theme_minimal()

# 4.3 UMAP-Plot
ggplot(
  umap_df,
  aes(
    UMAP1,
    UMAP2,
    label = sample,
    color = dex,
    shape = cell
  )
) +
  geom_point(size = 3) +
  geom_text(
    vjust = -0.7,
    size = 3
  ) +
  labs(
    title = "UMAP-Projektion",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Behandlung",
    shape = "Zelllinie"
  ) +
  theme_minimal()


# ============================================================
# Clustering
# ============================================================

# 1. Expressionsmatrix und Metadaten vorbereiten

metadata <- Probeninformationen

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


# 2. Euklidische Distanz und hierarchisches Clustering

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


# 3. Robustheitsprüfung verschiedener Linkage-Verfahren

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


# 4. Feature-Auswahl anhand der Genvarianz

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


# 5. Sensitivitätsanalyse der Feature-Anzahl

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


# 6. Sample-Distance-Heatmap

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


# 7. k-Means als unabhängiger Clustering-Ansatz

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


# 8. Silhouettenanalyse des k-Means-Ergebnisses

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
