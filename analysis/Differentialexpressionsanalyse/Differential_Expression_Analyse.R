# 1. Notwendigen Pakete für die Differential-Expression-Analyse laden
library (DESeq2)
library (tidyverse)
library (apeglm)

# 2. Daten aus csv laden
Daten_roh_long_filtered <- read_csv("data/Daten_roh_long_filtered.csv")
Tab_Probeninformationen <- read_csv("data/Tab_Probeninformationen.csv")
Tab_Probeninformationen <- as.data.frame(Tab_Probeninformationen)

# 3. Daten_roh_long_filtered für DESeq2 vorbereiten. Umwandlung in eine wide-format count Matrix. Die Spaltennamen
# der count matrix müssen mit den Zeilennamen der Metadaten übereinstimmen
count_mat <- Daten_roh_long_filtered %>%
  select (gene_id, sample, counts) %>%
  pivot_wider (
    names_from = sample,
    values_from = counts
  ) %>%
  column_to_rownames("gene_id") %>%
  as.matrix()
  

# 4. cell und dex als Faktoren definieren. Somit kann für die spätere Auswertung mit DESeq2 in Kategorien
# wie treated und untreated und in die vier unterschiedlichen Zelllinien/Spender eingeteilt werden.
Tab_Probeninformationen$cell <-
  factor(Tab_Probeninformationen$cell)

Tab_Probeninformationen$dex <-
  factor(Tab_Probeninformationen$dex)

rownames(Tab_Probeninformationen) <- Tab_Probeninformationen$Run

# Faktoren können überprüft werden
levels(Tab_Probeninformationen$cell)


# 5. Untreated als Referenzgruppe festlegen. So kann später intuitiver bestimmt werden, wie sich die Behandlung
# im Gegensatz zu keiner Behandlung auf die Genexpression auswirkt.
Tab_Probeninformationen$dex <- relevel(
  Tab_Probeninformationen$dex,
  ref = "untrt"
)

#Überprüfen ob untrt wirklich als Referenzlevel festgelegt wurde. Reihenfolge muss "untrt" "trt" sein
levels(Tab_Probeninformationen$dex)

# 6. Für DESeq2 müssen Reihenfolge der Count Matrix mit den Metadaten identisch sein. Dies wird überprüft. Falls es nicht stimmt
# bricht das Skript ab
stopifnot(
  identical(
    colnames(count_mat),
    rownames(Tab_Probeninformationen)
  )
)

# 7. DESeq2 durchführen
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

#8. Differentiale Genexpression zwischen behandelten (trt) und unbehandelten
# Proben (untrt) testen. Alpha = 0.05 beschreibt, dass Gene mit einem adjustierten p-value <0,05 signifikant sind.
res <- results(
  dds,
  contrast = c("dex","trt","untrt"),
  alpha = 0.05
)
summary(res)


# 9. Log2-FC shrinken, um unsichere oder möglicherweise überschätzte Effekte, besonders bei Genen mit wenigen Counts oder hoher Streuung, zu stabilisieren.
res_shrunk <- lfcShrink(
  dds,
  coef = "dex_trt_vs_untrt",
  type = "apeglm"
)

# 10. Ergebnisse sollen übersichtlicher gestaltet werden, weshalb Geninformationen dazugeholt werden. Zunächst wird eine Annotationstabelle erstellt.
annotation <- Daten_roh_long_filtered %>%
  select(
    gene_id,
    gene_name,
    symbol,
    gene_biotype
  ) %>%
  distinct()

# 11. DESeq2 Ergebnisse werden in Data Frame umgewandelt.
res_df <- as.data.frame(res_shrunk) %>%
  rownames_to_column("gene_id")

# 12. Annotation anhängen
res_df <- res_df %>%
  left_join(
    annotation,
    by = "gene_id"
  )

#13. Nach Signifikanz sortieren. Kleinsten padj kommen nach oben.
res_df <- res_df %>%
  arrange(padj)

#14. Eigene Tabelle für signifikante Gene (padj < 0.05) erstellen.
sig <- res_df %>%
  filter(
    !is.na(padj),
    padj < 0.05
  )

#15. Übersichtlichere Tabelle mit signifikanten Werten und wichtigen Informationen.
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

# 16. Betrachtung nur hochregulierter und dann nur herunterregulierter Gene nach der Behandlung.
sig_up <- sig %>%
  filter(log2FoldChange > 0)

sig_down <- sig %>%
  filter(log2FoldChange < 0)

# 17. Sortieren nach dem stärksten Effekt. Welche Gene besitzen einen hohen LFC. Welche wurden besonders stark hochreguliert?
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

#18. Welche Gene wurden am stärksten runterreguliert?
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

# 19.Ergebnistabelle als csv speichern
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

#20. MA-Plot aus geshrinkten LFCs
plotMA(res_shrunk, ylim=c(-2,2))

#Beim klicken auf einzelne Punkte im Plot können die zugehörigen Gene mit wichtigen Informationen angezeigt werden
idx <- identify(res_shrunk$baseMean, res_shrunk$log2FoldChange)
gene_ids <- rownames(res_shrunk)[idx]

res_df %>%
  filter(gene_id %in% gene_ids) %>%
  select(
    gene_id,
    symbol,
    gene_name,
    baseMean,
    log2FoldChange,
    padj
  )
