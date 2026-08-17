
# Allgemein

Im folgenden wird zunächst mittels dem Quality Assesment (QA) die Datenqualität überprüft. Anschließend werden mittels dem Quality Control spezifische Maßnahmen egriffen um die Datenqualität, basierend auf dem QA, zu verbessern. 

# QA

Bei dem ausgewählten Datensatz handelt es sich um RNA-Seq Daten. 
Zu Beginn wird die Library size überprüft. <img width="686" height="443" alt="Library size-final" src="https://github.com/user-attachments/assets/625570d4-fa5a-42a2-bf7c-ea6ae6f768a5" />
Dabei wird die Menge der counts in Abhängigkeit der Probe dargestellt. In diesem Balkendiagramm wird deutlich, dass SRR1039513 die geringste Größe aufweist. SRR1039517 dagegen hat die größte Library. Diese Erkenntnis deutet darauf hin, dass es unterschiedliche Sequenziertiefen gibt. Daher ist für zukünftige Auswertungen eine Normalisierung notwendig. Ebenfalls wird die Anzahl der detektierten Gene in Abhängigkeit der jeweiligen Probe dargestellt.
<img width="686" height="443" alt="anzahl gene-final" src="https://github.com/user-attachments/assets/348eb379-cd6f-4a36-81e8-3b37b6f0c3a5" />
Hier kann die selbe Beobachtung wie im Plot zuvor festgestellt werden. 

Danach wird eine Mean-Varianz-Analyse durchgeführt.
<img width="686" height="443" alt="Mean-varianz-plot-final" src="https://github.com/user-attachments/assets/d3012c26-7580-495f-951b-52bb33ba0fe0" />

Diese zeigt je höher die durchschnittliche Expression eines Gens ist, desto größer wird auch seine Varianz. Der Mean-Varianz-Plot der Rohdaten zeigt eine ausgeprägte Abhängigkeit der Varianz vom mittleren Expressionsniveau. Mit zunehmender mittlerer Countzahl steigt die Varianz der Gene deutlich an. Die Varianz liegt dabei deutlich über der erwarteten Varianz einer Poisson-Verteilung, was auf eine Overdispersion der RNA-seq-Countdaten hindeutet. Die Rohdaten weisen somit eine ausgeprägte heteroskedastische Struktur auf, weshalb für nachfolgende explorative Analysen eine Varianzstabilisierung sinnvoll ist. Anschließend wird ein Mean-Sd-Plot durchgeführt, welcher die Erkenntnis aus dem Mean-varianz-Plot bestätigt. 
<img width="686" height="441" alt="mean-sd-plot" src="https://github.com/user-attachments/assets/e7810a43-4cbc-4205-8440-770d6d1302c9" />
Die Standardabweichung hängt stark von der Expressionshöhe ab. Es wird anschließend eine log2-Transformation durchgeführt, weil die log2-Transformation die Abhängigkeit der Varianz vom Mittelwert reduziert. 
<img width="686" height="441" alt="mean-sd-log2" src="https://github.com/user-attachments/assets/0cf0171e-185b-49fd-8f0e-8da49b02ad3b" />Dieser Schritt soll in diesem Fall erstmal nur als Beobachtung genutzt werden um potentielle positive Effekte für das Preprocessing zu gewinnen. 

# QC

count 
