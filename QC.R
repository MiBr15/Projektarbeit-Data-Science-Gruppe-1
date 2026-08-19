# 1. PAKETE

library(edgeR)
library(dplyr)
library(tidyr)
library(readr)


# 2. DATEN EINLESEN

data <- read_csv(
  "Daten_roh_long.csv",
  show_col_types = FALSE
)


# 3. AUSGABEVERZEICHNIS

output_dir <- "C:/Users/komfo/OneDrive/Desktop/MASTER/2. Semester/limma_ausgabe"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# 4. SPALTEN CHECKEN

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


# 5. DATENTYPEN KONTROLLIEREN

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


# 6. COUNT-MATRIX AUS DER CSV ERSTELLEN

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


# 7. edgeR DGEList ERSTELLEN

dge <- DGEList(
  counts = count_matrix
)


# 8. LOW-COUNT-FILTERUNG MIT filterByExpr()

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


# 9. GEFILTERTE GENE WIEDER IN DIE LONG-TABELLE ZURÜCKFÜHREN

filtered_gene_ids <- rownames(
  filtered_count_matrix
)

data_filtered <- data %>%
  filter(
    gene_id %in% filtered_gene_ids
  )


# 10. GEFILTERTE DATEN ALS CSV SPEICHERN

output_file <- file.path(
  output_dir,
  "Daten_roh_long_filtered.csv"
)

write_csv(
  data_filtered,
  output_file,
  na = ""
)


# 11. AUSGABE FÜR ÜBERSICHT

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
