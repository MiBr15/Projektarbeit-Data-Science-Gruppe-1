
# Allgemein

Im Folgenden wird zunächst mithilfe des Quality Assessments (QA) die Qualität der Daten überprüft. Anschließend werden im Rahmen des Quality Controls (QC) spezifische Maßnahmen ergriffen, um die Datenqualität auf Grundlage der Ergebnisse des QA zu verbessern.
Bei dem ausgewählten Datensatz handelt es sich um RNA-Seq-Daten. Ziel dieses Aufgabenbereichs ist es, die Qualität des Datensatzes zu kontrollieren und anschließend gezielte Veränderungen am Datensatz vorzunehmen, um die Qualität zukünftiger Auswertungen zu verbessern. Dazu können beispielsweise zuvor identifizierte Ausreißer oder Gene mit sehr geringen Count-Werten entfernt werden.

Bei der Erstellung des Codes habe ich ChatGPT unterstützend genutzt. Insgesamt habe ich mich etwa zwei Tage intensiv mit der Projektarbeit beschäftigt und anschließend im weiteren Verlauf immer wieder kleinere Anpassungen und Verbesserungen vorgenommen.

# Quality Assesment
## Pakete installieren
Um im weiteren Verlauf die Daten analysieren und darstellen zu können müssen die Libraries "pheatmap" und "tidyverse" installiert werden. Somit können die Daten verarbeitet und dargestellt werden. 

## Datenstruktur kontrollieren
Zunächst wird die Datenstruktur, welche im Wrangling entstanden ist, nochmal auf ihren Aufbau untersucht. Dazu werden die Spalten-/Zeilennamen betrachtet, das Vorhandensein von Duplikaten oder fehlenden Werten und die menge an Daten die vorliegen. 

## Count-Matrix erstellen
Die Daten liegen im long Format vor. Für die Analyse der RNA-seq Daten wird allerdings eine Matrix aus Gene x Sample benötigt. Dieser Code-Abschnitt dreht somit einmal die Tabelle um, wodurch die Matrix entsteht. 

## Library size 
Zu Beginn wird die Library size überprüft. Dazu wird für jedes Sample der Gen-Count addiert. 
<img width="686" height="443" alt="Library size-final" src="https://github.com/user-attachments/assets/625570d4-fa5a-42a2-bf7c-ea6ae6f768a5" />

In diesem Balkendiagramm wird deutlich, dass SRR1039513 die geringste Library Größe aufweist. SRR1039517 dagegen hat die größte Library. Diese Erkenntnis deutet darauf hin, dass es unterschiedliche Sequenziertiefen gibt. Daher ist für zukünftige Auswertungen eine Normalisierung notwendig.

## Anzahl der detektierten Gene
Ebenfalls wird die Anzahl der detektierten Gene in Abhängigkeit der jeweiligen Probe dargestellt. Dazu wird spezifisch gezählz bei wie vielen Genen der Count größer als 0 ist.

<img width="686" height="443" alt="anzahl gene-final" src="https://github.com/user-attachments/assets/348eb379-cd6f-4a36-81e8-3b37b6f0c3a5" /> 

Hier kann die selbe Beobachtung wie im Plot zuvor festgestellt werden, dass SRR1039513 die geringste Anzahl an Genen aufweist. 

## Mean-Varianz-Analyse
Danach wird eine Mean-Varianz-Analyse durchgeführt.

<img width="686" height="443" alt="Mean-varianz-plot-final" src="https://github.com/user-attachments/assets/d3012c26-7580-495f-951b-52bb33ba0fe0" />

Diese zeigt je höher die durchschnittliche Expression eines Gens ist, desto größer wird auch seine Varianz. Der Mean-Varianz-Plot der Rohdaten zeigt eine ausgeprägte Abhängigkeit der Varianz vom mittleren Expressionsniveau. Mit zunehmender mittlerer Countzahl steigt die Varianz der Gene deutlich an. Die Varianz liegt dabei deutlich über der erwarteten Varianz einer Poisson-Verteilung, was auf eine Overdispersion der RNA-seq-Countdaten hindeutet. Die Rohdaten weisen somit eine ausgeprägte heteroskedastische Struktur auf, weshalb für nachfolgende explorative Analysen eine Varianzstabilisierung sinnvoll ist. 

## Mean-Sd-Plot
Anschließend wird ein Mean-Sd-Plot durchgeführt, welcher die Erkenntnis aus dem Mean-varianz-Plot bestätigt. 

<img width="686" height="441" alt="mean-sd-plot" src="https://github.com/user-attachments/assets/e7810a43-4cbc-4205-8440-770d6d1302c9" />

Die Standardabweichung hängt stark von der Expressionshöhe ab. 

## Log2 Transformation
Es wird anschließend eine log2-Transformation durchgeführt, weil die log2-Transformation die Abhängigkeit der Varianz vom Mittelwert reduziert. 

<img width="686" height="441" alt="mean-sd-log2" src="https://github.com/user-attachments/assets/0cf0171e-185b-49fd-8f0e-8da49b02ad3b" />

Dieser Schritt soll in diesem Fall erstmal nur als Beobachtung genutzt werden um potentielle positive Effekte für das Preprocessing zu gewinnen. 

# QC
Für das Quality Control beziehe ich mich auf die Erkenntnisse des QAs. Auffällig ist, dass im Datensatz einige Proben vorliegen mit einen count von Null oder einem sehr geringen Count. Somit wurde in diesem Schritt eine Filterung von low.counts und zero-counts durchgeführt. 

## Pakete laden
Für das QC werden die libraries "edgeR" verwendet  für die Filterung der Gene. "dplyr" kann die Daten direkt bearbeiten und mittels "tidyr" umformen. Zuletzt sollen die neu gefilterten daten in einer csv gespeicherrt werden, wodurch "readr" verwendet wird. 

## Daten kontrollieren
Zu Beginn wurde hierzu erneut die Datenstruktur kontrolliert. Da für das Eingreifen in die Daten wieder die Count-matrix benötigt wird, beziehe ich mich, auf die im QA gebildetete Matrix. 

## edgeR Daten erzeugen
Mittels DGEList können die Count-Daten und zusätzliche Informationen für zukünftige Analysen gespeichert werden. 

## Filter
Durch *filterByExpr(dge)* wird automatisch für jedes Gen analysiert ob es ausreichend stark exprimiert ist, um es für weitere Analysen im Datensatz zu behalten.  

## Gefilterte Tabelle
Nachdem zero-counts und low-counts entfernt wurden werden sie in die ursprüngliche Tabelle zurückgeführt und diese wird in einer csv Datei gespeichert. 




