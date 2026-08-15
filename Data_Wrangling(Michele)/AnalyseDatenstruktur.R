#1. Bioconductor, Paket und Datensatz laden

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("airway")

library(airway)

data(airway)

#2. Datenstruktur verstehen:

airway

assay(airway)

colData(airway)

head(assay(airway))

dim(airway)

#3. Datenextraktion

Messwerte_roh <- assay(airway)

Probeninformationen <- as.data.frame(colData(airway))

Geninformationen <- as.data.frame(rowData(airway))

#4. Untersuchung Messwerte_roh,Probeninformationen, Geninformationen

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

#5. Überprüfung gleicher Sample-IDs

all(colnames(Messwerte_roh) == Probeninformationen$Run)

#6. Erstellen einer Probeninformationstabelle

Tab_Probeninformationen <- Probeninformationen[, c(
  "Run",
  "cell",
  "dex",
  "avgLength"
)]

#7. Überprüfung erstellter Tabelle

nrow(Tab_Probeninformationen)

length(unique(Tab_Probeninformationen$Run))

#8. Messwerte_roh ins Long Format bringen

install.packages("tidyr")

library(tidyr)

Messwerte_roh_long <- as.data.frame(Messwerte_roh)

Messwerte_roh_long$gene_id <- rownames(Messwerte_roh_long)

rownames(Messwerte_roh_long) <- NULL

Messwerte_roh_long <- pivot_longer(
  Messwerte_roh_long,
  col = -gene_id,
  names_to ="sample",
  values_to = "counts"
  
)

#9. Überprüfung der Tabelle

dim(Messwerte_roh_long)

#10. Messwerte_roh_long mit Tab_Probeninformationen verbinden

library(dplyr)

Daten_roh_long <- Messwerte_roh_long %>%
  left_join(
    Tab_Probeninformationen,
    by = c("sample" = "Run")
  )

#11. Neue Tabelle überprüfen

head(Daten_roh_long)

dim(Daten_roh_long)

sum(is.na(Daten_roh_long$cell))


#12. Join von Daten_roh_long und Geninformationen

Daten_roh_long <- Daten_roh_long %>%
  left_join(
    Geninformationen,
    by = "gene_id"
  )

#13. Überprüfung des Daten_roh_long-Datensatzes

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

#14. Rohdatentabelle als csv speichern

write.csv(
  Daten_roh_long,
  "Daten_roh_long.csv",
  row.names = FALSE
)
