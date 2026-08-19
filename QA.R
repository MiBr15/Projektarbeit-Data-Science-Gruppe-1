# 1. PAKETE
#install.packages("tidyverse")
#install.packages("pheatmap")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

library(tidyverse)
library(pheatmap)


# 2. csv einlesen 

Daten_roh_long <- read.csv(
  "Daten_roh_long.csv",
  stringsAsFactors = FALSE
)

# Welche Daten verwendet werden
data <- Daten_roh_long


# 3. KONTROLLE DER DATENSTRUKTUR

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


# 4. Count-Matrix

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


# 5. Library size

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


# 6. Anzahl der detektierten Gene

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


# 7. Mean-Varianz-Analyse

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


# 8. Mean-SD-Plot

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


# 9. Log2 Transformation und Mean-SD-Plot mit log2 Daten

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
