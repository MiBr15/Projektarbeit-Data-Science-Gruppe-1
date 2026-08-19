library(edgeR)
library(dplyr)
library(tidyr)
library(readr)

# Daten
data <- Daten_roh_long

output_dir <- "C:/Users/komfo/OneDrive/Desktop/MASTER/2. Semester/limma_ausgabe"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Spalten checken
required_cols <- c("gene_id", "sample", "counts")

missing_cols <- setdiff(required_cols, colnames(data))

if (length(missing_cols) > 0) {
  stop(
    "Fehlende Pflichtspalten: ",
    paste(missing_cols, collapse = ", ")
  )
}

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


count_matrix 

# 1. edgeR DGEList erstellen

dge <- DGEList(counts = count_matrix)


# 2. Low-Count-Filterung mit filterByExpr()

keep <- filterByExpr(dge)

cat("\nGene vor Filterung: ", nrow(dge), "\n")
cat("Gene nach Filterung: ", sum(keep), "\n")
cat("Entfernte Gene: ", sum(!keep), "\n")

if (sum(keep) == 0) {
  stop("Kein Gen erfüllt die filterByExpr()-Kriterien.")
}

filtered_count_matrix <- count_matrix[keep, , drop = FALSE]


# 3. Gefilterte Gene wieder in die ursprüngliche Long-Tabelle zurückführen

filtered_gene_ids <- rownames(filtered_count_matrix)

data_filtered <- data %>%
  filter(gene_id %in% filtered_gene_ids)

# 4. als csv speichern 
timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

output_file <- file.path(
  output_dir,
  paste0("Daten_roh_long_filtered_", timestamp, ".csv")
)

write_csv(
  data_filtered,
  output_file,
  na = ""
)

# Ausgabe für Übersicht

cat("\nLow-Count-Filterung abgeschlossen.\n")
cat("Gene vor Filterung: ", nrow(count_matrix), "\n")
cat("Gene nach Filterung: ", nrow(filtered_count_matrix), "\n")
cat("Entfernte Gene: ", sum(!keep), "\n")
cat("Zeilen vor Filterung: ", nrow(data), "\n")
cat("Zeilen nach Filterung: ", nrow(data_filtered), "\n")
cat("Ausgabedatei: ", output_file, "\n")