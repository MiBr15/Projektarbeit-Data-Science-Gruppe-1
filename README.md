# Fragestellung
Mithilfe der differential expression analysis kann untersucht werden, wie sich die Genexpression zwischen unterschiedlichen experimentellen Gruppen verändert. Im Beispiel des airway Datensatzes kann mithilfe dieser Analyse ein Vergleich der Genexpression zwischen unbehandelten Zellen und mit Dexamethason behandelten Zellen stattfinden. Dabei wird für jedes Gen unter anderem ein log2FoldChange (log2FC) berechnet, der die Stärke und Richtung der Expressionsänderung beschreibt. Ein positiver log2FC bedeutet eine höhere Expression unter der Behandlung, ein negativer log2FC eine niedrigere Expression. Zusätzlich kann bestimmt werden, ob diese Unterschiede statistisch signifikant sind.
# Vorbereitung der Daten
Für die Differential-Expression-Analyse wird das DESeq2 package verwendet. Hierfür müssen die Daten nach dem data wrangling und dem QC zunächst noch weiter vorbereitet werden, bevor das package die Analyse durchführen kann. Die gefilterten Daten dürfen nicht normalisiert vorliegen, da DESeq2 eine eigene Normalisierung durchführt. Die Ausgangsdaten liegen im Long-Format vor, sodass jedes Gen für jedes Sample eine eigene Zeile besitzt. Mit pivot_wider() werden die einzelnen Samples deshalb zu Spalten umgeformt. Die Spaltennamen der Count-Matrix müssen außerdem mit den Zeilennamen der Metadaten übereinstimmen.
# Faktoren bestimmen
Damit DESeq2 die experimentellen Gruppen unterscheiden kann, werden dex und cell als Faktoren definiert. Der Faktor dex unterscheidet zwischen unbehandelten (untrt) und behandelten (trt) Proben, während cell die vier unterschiedlichen Zelllinien bzw. Spender beschreibt. Anschließend wird untrt als Referenzlevel für dex festgelegt. Dadurch beschreibt ein positiver log2FoldChange später eine höhere Expression unter Dexamethason-Behandlung im Vergleich zu unbehandelten Zellen. 
# DESeq2 ausführen
Kurz bevor DESeq2 ausgeführt wird, wird nochmals überprüft, ob die Spaltennamen der Count-Matrix mit den Zeilennamen der Metadaten übereinstimmen, da dies essenziell für DESeq2 ist. Wenn dies nicht der Fall sein sollte, bricht das Skript automatisch ab, bevor DESeq2 ausgeführt wird. Anschließend wird mit den Informationen der Count-Matrix und der Metadaten DESeq2 ausgeführt. Dadurch werden Unterschiede zwischen den Zelllinien/Spendern im statistischen Modell berücksichtigt, sodass der Effekt der Dexamethason-Behandlung möglichst unabhängig davon geschätzt werden kann. Bei der DESeq2 laufen folgende Schritte ab:
- Normalisierung der unterschiedlichen Sequenziertiefen
- Schätzung der Dispersion
- Anpassung des statistischen Modells
- Test auf differentielle Genexpression
Zuletzt wird überprüft, ob tatsächlich trt zu utrt relativiert wurde. Wenn dies nicht der Fall sein sollte, bricht das Skript automatisch ab.
# Vergleich zwischen behandelten und nicht behandelten Zellen
Anschließend werden die Ergebnisse für den Vergleich trt gegenüber untrt extrahiert. Für jedes Gen werden unter anderem der baseMean als mittlere normalisierte Count-Zahl über alle Samples, der log2FoldChange als Stärke und Richtung der Expressionsänderung, der p-Wert sowie der für multiples Testen korrigierte p-Wert (padj) ausgegeben. Gene mit padj < 0.05 werden in dieser Analyse als statistisch signifikant betrachtet.
# log2-FC stabilisieren
Da der log2FoldChange insbesondere bei Genen mit geringer Information oder hoher Streuung unsicher und teilweise überschätzt sein kann, werden die log2FoldChanges mit lfcShrink() stabilisiert. Unsichere Schätzungen werden dabei stärker in Richtung null gezogen, während gut durch die Daten gestützte Effekte meist nur gering verändert werden. Das Shrinkage verbessert dadurch insbesondere die Interpretation, Darstellung und das Ranking der Effektgrößen. Die statistische Signifikanz (pvalue und padj) wird dadurch nicht verändert.
# Aufarbeitung der Darstellung
Um schlussendlich die Ergebnisse übersichtlich und leicht interpretierbar darzustellen, wird zunächst eine Annotationstabelle erstellt, die wichtige Informationen zu den einzelnen Genen enthält. Um diese Annotationen an die Ergebnisse anzufügen, müssen diese Ergebnisse zunächst mithilfe von tidyverse von einem Bioconductor Objekt in einen normalen Data Frame umgewandelt werden. Anschließend können beide Tabellen mit leftjoin() über die identischen gene ids zusammengefügt werden. Das Ergebnis ist eine übersichtliche Tabelle mit den statistischen Daten und den zugehörigen Informationen der Gene. 
In den darauffolgenden Schritten kann diese Tabelle je nach Verwendungszweck sortiert oder einzelne Werte extrahiert werden. Mit dem Skript können folgende Tabellen erstellt werden:
- sortiert nach höchster Signifikanz
- nur signifikante Werten (padj<0.05)
- nur die wichtigsten Gen- und statistischen Informationen
- nur hoch- oder nur runterregulierte Gene
- die am stärksten hoch- und die am stärksten runterregulierten Gene
Dabei stellt sich heraus, dass nach der Behandlung das Gen SPARCL1 mit der gene id ENSG00000152583 einenlog2-FC von 4,55 besitzt und am stärksten hochreguliert wurde. Dahingegen besitzt das Gen VCAM1 mit der gene id ENSG00000162692 einen log2-FC von -3,68 und wurde am stärksten runterreguliert. 
Die gesamte annotierte Tabelle und die Tabelle mit nur den signifikanten Ergebnisse werden als besonders wichtig betrachtet, sodass diese als .csv gespeichert werden können.
# Interaktiver MA-plot
Anschließend wird ein MA-Plot aus den geshrinkten log2FoldChanges erstellt. Für jedes Gen wird auf der x-Achse der baseMean, also die mittlere normalisierte Count-Zahl, und auf der y-Achse der geshrinkte log2FoldChange dargestellt. Gene oberhalb von 0 sind unter Dexamethason höher exprimiert, Gene unterhalb von 0 niedriger exprimiert. 
<img width="1042" height="884" alt="MA-plot" src="https://github.com/user-attachments/assets/1186ce5e-06ec-413b-a468-45d9818f817b" />
Mithilfe von identify() können einzelne Punkte anschließend interaktiv ausgewählt und die zugehörigen Gene zusammen mit ihren biologischen und statistischen Informationen angezeigt werden.    
