#1. Laden aller  nötigen Pakete zum Preprocessing

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2", ask = FALSE, update = FALSE)
}

library(DESeq2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

#2. Einlesen der gefilterten Rohdaten

fil_Daten <- read.csv("Daten_roh_long_filtered.csv")

#3. Prüfung der überreichten Daten_roh_long_filtered.csv

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

#4. Extrahierung einer Count-Matrix aus fil_Daten

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

#5. Probeninformationen (Metadaten) zuordnen

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

#6. Beginn des Preprocessing/Normalisierens - Erstellen eines DESeq2-Objekts

Probeninformationen$cell <- factor(Probeninformationen$cell)
Probeninformationen$dex <- factor(Probeninformationen$dex)

dds <- DESeqDataSetFromMatrix(
  countData = Daten_matrix,
  colData = Probeninformationen,
  design = ~ cell + dex
)

#7. Size-Factor-Normalisierung

dds <- estimateSizeFactors(dds)

sizeFactors(dds)

#8. Normalisierung der Messdaten

normalisierte_counts <- counts(
  dds,
  normalized = TRUE
)

#9. Kontrolle der normalisierten Counts

dim(normalisierte_counts)

head(normalisierte_counts[, 1:4])

data.frame(
  sample = colnames(Daten_matrix),
  size_factor = sizeFactors(dds)
)

#10. Variance Stabilizing Transformation (VST)

vst_daten <- vst(
  dds,
  blind = TRUE
)

vst_matrix <- assay(vst_daten)

#11. QA der Size Factors

size_factor_info <- data.frame(
  sample = colnames(Daten_matrix),
  cell = Probeninformationen$cell,
  dex = Probeninformationen$dex,
  size_factor = sizeFactors(dds)
)

print(size_factor_info)

#12. QA der Library Size

library_größe <- colSums(Daten_matrix)

library_info <- data.frame(
  sample = names(library_größe),
  library_größe = as.numeric(library_größe),
  size_factor = sizeFactors(dds)
)

print(library_info)

#13. QA der VST-Matrix

dim(vst_matrix)

head(vst_matrix[, 1:4])

#14. QA des Mean-Variance-Plots

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

#15. Ergebnisse speichern

output_dir <- getwd()

vst_ergebnisse <- as.data.frame(vst_matrix) %>%
  rownames_to_column("gene_id")

normalisierte_ergebnisse <- as.data.frame(normalisierte_counts) %>%
  rownames_to_column("gene_id")

normalisierungs_info <- Probeninformationen %>%
  mutate(
    size_factor = sizeFactors(dds)
  )

png(
  filename = file.path(output_dir, "VST_Mittelwert-Varianz-Beziehung.png"),
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
    output_dir,
    paste0("DESeq2_normalisierte_counts.csv")
  ),
  na = ""
)

write_csv(
  vst_ergebnisse,
  file.path(
    output_dir,
    paste0("VST_transformierte_Daten.csv")
  ),
  na = ""
)

write_csv(
  normalisierungs_info,
  file.path(
    output_dir,
    paste0("Probeninformationen_Normalisierung.csv")
  ),
  na = ""
)

message("Alles fertig")
