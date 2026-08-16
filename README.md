# Allgemein

Im folgenden wird zunächst mittels dem Quality Assesment (QA) die Datenqualität überprüft. Anschließend werden mittels dem Quality Control spezifische Maßnahmen egriffen um die Datenqualität, basierend auf dem QA, zu verbessern. 

# QA

Bei dem ausgewählten Datensatz handelt es sich um RNA-Seq Daten. Daher wird zunächst die Verteilung der Counts betrachtet. Dabei fällt auf, dass die acht Samples sich relativ ähnlich verhalten. Die Daten wiesen erwartungsgemäß eine stark rechtsschiefe Verteilung mit einem hohen Anteil niedrig exprimierter Gene und wenigen Genen mit hohen Counts auf.
Anschließend die Count Verteilung pro sample. Hier ist auffällig, dass die Probe SRR1039513 im Vergleich zu den anderen Proben weniger detektierte Gene aufweist. Diese Daten werden dann als eine Matrix dargestellt. 
Die Library size wird ebenfalls überprüft. In diesem Balkendiagramm wird deutlich, dass SRR1039513 die geringste Größe aufweist. SRR1039517 dagegen hat die größte Library. Diese Erkenntnis deutet darauf hin, dass es unterschiedliche Sequenziertiefen gibt. Daher ist für weitere Auswertungen eine Normalisierung notwendig.  

Danach wird eine Mean-Varianz-Analyse durchgeführt. Diese zeigt je höher die durchschnittliche Expression eines Gens ist, desto größer wird auch seine Varianz. Der Mean-Varianz-Plot der Rohdaten zeigt eine ausgeprägte Abhängigkeit der Varianz vom mittleren Expressionsniveau. Mit zunehmender mittlerer Countzahl steigt die Varianz der Gene deutlich an. Die Varianz liegt dabei deutlich über der erwarteten Varianz einer Poisson-Verteilung, was auf eine Overdispersion der RNA-seq-Countdaten hindeutet. Die Rohdaten weisen somit eine ausgeprägte heteroskedastische Struktur auf, weshalb für nachfolgende explorative Analysen eine Varianzstabilisierung sinnvoll ist. Anschließend wird ein Mean-Sd-Plot durchgeführt, welcher die Erkenntnis aus dem Mean-varianz-Plot bestätigt. Die Standardabweichung hängt stark von der Expressionshöhe ab. Es wird anschließend eine log2-Transformation durchgeführt, weil die log2-Transformation die Abhängigkeit der Varianz vom Mittelwert reduziert.

# QC

Heatmap auch sinvoll um zwischen den Samples potentielle Ausreißer zu identifizierne. 

Die PCA auf Basis der VST-transformierten RNA-seq-Daten zeigt eine deutliche Separation der behandelten und unbehandelten Samples entlang der ersten Hauptkomponente (PC1), welche 47,7 % der Gesamtvarianz erklärt. Die ersten beiden Hauptkomponenten erklären zusammen 71,2 % der Gesamtvarianz. Die klare Trennung der Samples nach Behandlungsstatus deutet darauf hin, dass die Behandlung einen wesentlichen Einfluss auf die globalen Genexpressionsprofile hat. Gleichzeitig zeigen sich Unterschiede zwischen den Zelllinien, insbesondere entlang der zweiten Hauptkomponente (PC2, 23,5 %). Ein offensichtlicher Sample-Outlier ist im PCA nicht erkennbar. Die Replikate weisen innerhalb der jeweiligen experimentellen Gruppen ein insgesamt konsistentes Muster auf.
