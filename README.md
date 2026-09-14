# ReadMe - Preprocessing_Code.R

## Erklärung zur KI-Nutzung

Zur Erstellung des Codes „Preprocessing_Code.R“ wurde die ChatGPT-Version „GPT-5.6 Luna“ unterstützend eingesetzt. 
Die Nutzung beschränkte sich dabei auf:

- Die Erstellung und Korrektur der Codesyntax
- Umsetzung und Korrektur einzelner Codeabschnitte
- Klärung von Fragen zu R-Funktionen und Paketen

Der grundlegende Aufbau sowie die verwendeten Konzepte des Codes wurden vorgegeben. Sämtliche von der KI vorgeschlagenen Codeabschnitte wurden auf fachliche Richtigkeit sowie ihre Sinnhaftigkeit überprüft und nachvollzogen.

## Allgemein

Der in diesem Branch vorliegende Code "Preprocessing_Code.R" dient der weiteren Prozessierung der Rohdaten nach erfolgter Low-Count-Filterung im Branch "QC-Emily". Hierbei werden die gefilterten Rohdaten, nach einer Kontrolle der Datenqualität und -integrität, mittels des Bioconductor-Pakets DESeq2 einer Size-Factor-Normalisierung unterzogen und nachfolgend zudem eine Variance Stabilizing Transformation durchgeführt.

Der folgende Abschnitt enthält eine detaillierte Erklärung zur Funktionsweise des für diese Aufgabe erstellten Codes. 

## Struktur und Funktion des Codes

### 1. Laden aller nötigen Pakete zum Preprocessing

In diesem Abschnitt wird, sofern nötig, der BiocManager installiert und nachfolgend die notwendigen R-Pakete geladen.

- DESeq2 (Normalisierung und VST)
- dplyr (Datenmanipulation)
- tydyr (Umformung von Datenstrukturen)
- readr (Einlesen und Speichern von csv-Dateien)
- tibble (Verarbeitung der Data Frames und Row Names)

### 2. Einlesen der gefilterten Rohdaten

In diesem Schritt wird die Datei "Daten_roh_long_filtered.csv" eingelesen und im Objekt fil_Daten gespeichert. Hierbei handelt es sich um die gefilterten Rohdaten im Long-Format, die im vorausgehenden Schritt im Branch "QC-Emily" erzeugt worden sind.

### 3. Prüfung der überreichten Daten_roh_long_filtered.csv

Im Folgenden findet eine nochmalige Kontrolle der in der csv-Datei überreichten Daten statt. Hierbei werden folgende Aspekte untersucht.

- Prüfung der benötigten Spalten (gene_id, sample, Counts, cell, dex)
- Prüfung auf NA-Werte, negative Counts oder nicht-ganzzahlige Einträge
- Prüfung ob ob jede Probe (sample) eindeutig einem Zelltyp (cell) und einer Behandlung (dex) zugeordnet werden kann. 

Sollte eine dieser Prüfungen eine Fehlermeldung hervorbringen, wird die Analyse abgebrochen.

### 4. Extrahierung einer Count-Matrix aus fil_Daten

Zur Anwendung von DESeq2 werden die gefilterten Rohdaten nicht im Long-Format, sondern als Count-Matrix benötigt. Diese Count-Matrix (Daten_matrix) wird in diesem Codeabschnitt aus fil_Daten generiert. Die Gene werden als Row Names und Samples als Column Names genutzt. Die Matrix wird als Integer-Matrix gespeichert und die Datenstruktur nochmals überprüft.

### 5. Probeninformationen (Metadaten) zuordnen

Neben der Count-Matrix werden aus den gefilterten Rohdaten ebenfalls die Metadaten extrahiert und in "Probeninformationen" gespeichert. Die Metadaten werden nachfolgend in dieselbe Reihenfolge gebracht wie die Spalten der Count-Matrix. Am Ende wird geprüft ob jedes Sample den richtigen Metadaten zugeordnet werden kann. Sollte dies nicht möglich sein, wird die Analyse gestoppt

### 6. Beginn des Preprocessing/Normalisierens - Erstellen eines DESeq2-Objekts

Zu Beginn des Preprocessing wird in diesem Codesegment das DESeq2-Objekt erstellt. Hierfür werden "cell" und "dex" als Faktoren definiert und anschließend im DESeq2-Objekt "DESeqDataSet die Rohcounts, die experimentellen Metadaten und das experimentelle Design miteinander verknüpft.

### 7. Size-Factor-Normalisierung

In diesem Abschnitt wird die Library-Size-Normalisierung von DESeq2 durchgeführt. Mittels "estimateSizeFactors(dds)" wird für jedes Sample unter Berücksichtigung der Sequenziertiefe und Library Size ein Size Factor berechnet und anschließend durch "sizeFactors(dds)" ausgegeben

### 8. Normalisierung der Messdaten

Mit den vorher bestimmten Size Factors werden in diesem Abschnitt die Count-Werte um Unterschiede in der Sequenziertiefe und Library normalisiert und erneut als Matrix "DESeq2_normalisierte_counts" (gleiches Format wie die Count-Matrix) ausgegeben.

### 9. Kontrolle der normalisierten Counts

Die Dimension, die ersten Werte und die berechneten Size Factors werden hier nochmals auf Vollständigkeit der Gene und Samples sowie auf stattgefundene Normalisierung kontrolliert.

### 10. Variance Stabilizing Transformation (VST)

In diesem Abschnitt wird eine Variance Stabilizing Transformation (VST) durchgeführt um heteroskedastische Muster bei den Daten zu entfernen und diese somit für statische und graphische Analysen vorzubereiten. Die so vorbereiteten Daten sind als "vst_matrix" gespeichert.

### 11. QA der Size Factors

Zur Überprüfung der Size Factors werden Sample, Zelltyp, Behandlung und Size Factor zusammengeführt und visuell geprüft ob bestimmte Samples auffällig hohe oder niedrige Size Factors besitzen. Stark abweichende Size Factors können stark unterschiedliche Library Sizes hinweisen und sollten nochmals geprüft werden. 

### 12. QA der Library Size

In diesem Abschnitt wird die ursprüngliche Library Size jedes Samples bestimmt und zusammen mit dem Size Factor und Sample ausgegeben. Durch Betrachtung dieser Informationen können Unterschiede zwischen der ursprünglichen Sequenziertiefe und den für die Normalisierung berechneten Size Factors näher analysiert werden.

### 13. QA der VST-Matrix

Zur Bestimmung einer erfolgreich durchgeführten VST-Transformation werden die Dimension und die ersten Werte der VST-Matrix ausgegeben und auf Richtigkeit überprüft.

### 14. QA des Mean-Variance-Plots

Dieser Abschnitt untersucht die Mean-Variance-Beziehung nach der VST. Hierfür wird zunächst die mittlere VST-Expression über alle Samples und die Varianz aller VST-Werte bestimmt und in einem Scatterplot dargestellt.

<img width="1600" height="1200" alt="VST_Mittelwert-Varianz-Beziehung" src="https://github.com/user-attachments/assets/ce8667b4-95c9-4635-872b-f15db44a46a5" />


In der Graphik ist deutlich zu erkennen, dass die in der QC festgestellte Heteroskedastizität bedeutend abgeschwächt wurde und die VST-Transformation damit erfolgreich war.

### 15. Ergebnisse speichern

Im finalen Schritt werden die erstellten Daten als csv-Dateien oder png-Graphiken gespeichert.

- VST_transformierte_Daten.csv (Enthält die VST-Matrix; zu nutzen für explorative Analysen, PCA, Clustering oder Heatmaps)

- DESeq2_normalisierte_counts.csv (Enthält die durch die Size Factors normalisierten Counts; in dieser Analyse nur als Zwischenergebnis erstellt)

- Probeninformationen_Normalisierung.csv (Enthält die Sample-Metadaten und die berechneten Size Factors; dient der Nachvollziehbarkeit der Normalisierung)

- VST_Mittelwert-Varianz-Beziehung.png (Enthält die Darstellung des in Punkt 14 durchgeführten QA der VST)
