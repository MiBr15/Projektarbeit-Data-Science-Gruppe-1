# Data-Science Projekt – Gruppe 1

## Projektübersicht

In diesem Projekt wird der RNA-Seq-Datensatz [`airway`](https://bioconductor.org/packages/airway/) aus Bioconductor analysiert. Der Datensatz umfasst acht Proben aus vier primären humanen Airway-Smooth-Muscle-Zelllinien. Für jede Zelllinie liegt jeweils eine unbehandelte Probe (`untrt`) und eine mit Dexamethason behandelte Probe (`trt`) vor.

Der Datensatz basiert auf der Publikation:

> Himes et al. (2014): *RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells*.  
> DOI: https://doi.org/10.1371/journal.pone.0099625

Ziel der Analyse ist es, die Daten reproduzierbar aufzubereiten, ihre Qualität zu untersuchen und zu prüfen, welche Unterschiede zwischen behandelten und unbehandelten Proben erkennbar sind.

## Analyseschritte

Die Auswertung umfasst:

1. **Data Wrangling**  
   Aufbereitung der Count-Daten sowie Verknüpfung mit Proben- und Geninformationen.

2. **Quality Assessment und Quality Control**  
   Untersuchung von Library Size, Anzahl detektierter Gene sowie Mean-Variance- und Mean-SD-Beziehungen. Anschließend werden Gene mit sehr geringer Expression gefiltert.

3. **Preprocessing und Normalisierung**  
   Normalisierung mit DESeq2 und Variance Stabilizing Transformation (VST).

4. **Differential-Expression-Analyse**  
   Untersuchung Dexamethason-assoziierter Expressionsänderungen mit einem DESeq2-Modell (`~ cell + dex`) und LFC-Shrinkage mit `apeglm`.

5. **Dimensionsreduktion und Visualisierung**  
   Explorative Darstellung der Proben mittels PCA, t-SNE und UMAP.

6. **Clustering**  
   Hierarchisches Clustering, Vergleich verschiedener Linkage-Methoden, Sensitivitätsanalyse bezüglich der verwendeten Gene sowie k-Means-Clustering und Silhouettenanalyse.

Aufgrund der kleinen Stichprobengröße von acht Proben werden insbesondere t-SNE und UMAP ausschließlich explorativ interpretiert.

## Projektstruktur

```text
.
├── README.md
├── report_Gruppe_1.qmd
├── data/
├── images/
└── analysis/
    ├── Clustering/
    ├── Data_Processing/
    ├── Data_Wrangling/
    ├── Differentialexpressionsanalyse/
    ├── Dimensionsreduktion_und_Visualisierung/
    └── QA-QC/
```

Die einzelnen Unterordner enthalten die jeweils zugehörigen Analyse-Skripte und ergänzenden README-Dateien. Die zentrale reproduzierbare Gesamtauswertung befindet sich in `report_Gruppe_1.qmd`.

## Abhängigkeiten

Für die Analyse werden insbesondere folgende R-Pakete verwendet:

### Bioconductor

- `airway`
- `DESeq2`
- `edgeR`
- `apeglm`

### CRAN

- `tidyverse`
- `pheatmap`
- `Rtsne`
- `uwot`
- `cluster`

Installation der Bioconductor-Pakete:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(
  c("airway", "DESeq2", "edgeR", "apeglm")
)
```

Installation der CRAN-Pakete:

```r
install.packages(
  c("tidyverse", "pheatmap", "Rtsne", "uwot", "cluster")
)
```

## Reproduzierbarkeit

Das Repository sollte als Arbeitsverzeichnis geöffnet werden, damit die relativen Pfade zu `data/` und `images/` funktionieren.

Der vollständige Bericht kann mit Quarto aus `report_Gruppe_1.qmd` erzeugt werden. Die R-Code-Chunks enthalten die vollständigen Analyseschritte; erläuternder Text und Interpretation stehen außerhalb der Code-Chunks.

## Autor:innen

- Michele Brömme
- Felix Kieborz
- Emily Komfort
- Felix Meyer
- Rebecca Schreck

## Nutzung von KI-Werkzeugen

Für einzelne Teile des Projekts wurden KI-Werkzeuge unterstützend eingesetzt, beispielsweise zur Erstellung und Korrektur von R-Syntax, zur Klärung von Fragen zu R-Funktionen und Paketen sowie zur Überprüfung einzelner Analyseschritte.

Die Auswahl der grundlegenden Analysekonzepte und Methoden erfolgte durch die Projektgruppe. KI-generierte Vorschläge wurden überprüft, angepasst und in die gemeinsame reproduzierbare Analyse integriert.
