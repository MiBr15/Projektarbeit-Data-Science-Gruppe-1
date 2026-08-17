# Allgemein

Im folgenden wird zunächst mittels dem Quality Assesment (QA) die Datenqualität überprüft. Anschließend werden mittels dem Quality Control spezifische Maßnahmen egriffen um die Datenqualität, basierend auf dem QA, zu verbessern. 

# QA

Bei dem ausgewählten Datensatz handelt es sich um RNA-Seq Daten. 
Zu Beginn wird die Library size überprüft. dabei wird die Menge der counts in Abhängigkeit der Probe dargestellt. In diesem Balkendiagramm wird deutlich, dass SRR1039513 die geringste Größe aufweist. SRR1039517 dagegen hat die größte Library. Diese Erkenntnis deutet darauf hin, dass es unterschiedliche Sequenziertiefen gibt. Daher ist für zukünftige Auswertungen eine Normalisierung notwendig. Ebenfalls wird die Anzahl der detektierten Gene in Abhängigkeit der jeweiligen Probe dargestellt. Hier kann die selbe Beobachtung wie im Plot zuvor festgestellt werden. 

Danach wird eine Mean-Varianz-Analyse durchgeführt. Diese zeigt je höher die durchschnittliche Expression eines Gens ist, desto größer wird auch seine Varianz. Der Mean-Varianz-Plot der Rohdaten zeigt eine ausgeprägte Abhängigkeit der Varianz vom mittleren Expressionsniveau. Mit zunehmender mittlerer Countzahl steigt die Varianz der Gene deutlich an. Die Varianz liegt dabei deutlich über der erwarteten Varianz einer Poisson-Verteilung, was auf eine Overdispersion der RNA-seq-Countdaten hindeutet. Die Rohdaten weisen somit eine ausgeprägte heteroskedastische Struktur auf, weshalb für nachfolgende explorative Analysen eine Varianzstabilisierung sinnvoll ist. Anschließend wird ein Mean-Sd-Plot durchgeführt, welcher die Erkenntnis aus dem Mean-varianz-Plot bestätigt. Die Standardabweichung hängt stark von der Expressionshöhe ab. Es wird anschließend eine log2-Transformation durchgeführt, weil die log2-Transformation die Abhängigkeit der Varianz vom Mittelwert reduziert. Dieser Schritt soll in diesem Fall erstmal nur als Beobachtung genutzt werden um potentielle positive Effekte für das Preprocessing zu gewinnen. 

# QC

count 
