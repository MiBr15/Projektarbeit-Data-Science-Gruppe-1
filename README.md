# Clustering

## Fragestellung und Hypothese

Mittels unüberwachtem Clustering wurde untersucht, ob sich die acht Airway-Proben allein anhand ihrer Genexpressionsprofile in biologisch sinnvolle Gruppen einteilen lassen. Dabei wurde insbesondere geprüft, ob die resultierende Clusterstruktur mit der Behandlung durch Dexamethason oder mit den vier untersuchten Zelllinien übereinstimmt.

Es wurde die Hypothese aufgestellt, dass die Behandlung mit Dexamethason zu systematischen Veränderungen des Genexpressionsprofils führt und sich behandelte und unbehandelte Proben daher auch ohne Verwendung der bekannten Gruppeninformation voneinander trennen lassen.

## Verwendete Daten

Für die Clustering-Analyse wurden die nach der Qualitätskontrolle verbliebenen und mittels Variance Stabilizing Transformation (VST) transformierten Expressionsdaten verwendet. Die Expressionsmatrix enthält 14.224 Gene für acht Proben.

Benötigte Dateien:

- `VST_transformierte_Daten.csv`
- `Probeninformationen_Normalisierung.csv`

## Hierarchisches Clustering

Zur Quantifizierung der Ähnlichkeit zwischen den Proben wurde die euklidische Distanz anhand der VST-transformierten Expressionswerte berechnet. Anschließend wurde ein hierarchisches Clustering mit Complete Linkage durchgeführt.

![Hierarchisches Clustering der Airway-Proben](images/clustering_complete.png)

Das hierarchische Clustering führte zu zwei klar getrennten Hauptclustern. Alle vier unbehandelten Proben (`untrt`) wurden einem Cluster und alle vier mit Dexamethason behandelten Proben (`trt`) dem zweiten Cluster zugeordnet. Die Information über den Behandlungsstatus wurde nicht für die Berechnung der Cluster verwendet.

Innerhalb beider Behandlungsgruppen zeigten insbesondere die Proben der Zelllinien N052611 und N061011 eine hohe Ähnlichkeit. Die Probe der Zelllinie N080611 wies innerhalb beider Behandlungsgruppen eine vergleichsweise größere Distanz zu den übrigen Proben auf.

## Einfluss des Linkage-Verfahrens

Zur Überprüfung der Robustheit des Ergebnisses wurde das hierarchische Clustering zusätzlich mit Average Linkage und Ward.D2 durchgeführt.

Bei allen drei untersuchten Linkage-Verfahren – Complete, Average und Ward.D2 – entstanden dieselben beiden Hauptgruppen entsprechend dem Behandlungsstatus. Die Trennung zwischen behandelten und unbehandelten Proben war damit gegenüber der Wahl des Linkage-Verfahrens robust.

## Einfluss der Feature-Auswahl

Zusätzlich wurde untersucht, ob die Clusterstruktur von der Anzahl der verwendeten Gene abhängt. Hierzu wurden die Gene anhand ihrer Varianz über die acht Proben geordnet. Die Auswahl erfolgte unabhängig vom bekannten Behandlungsstatus.

Zunächst wurden die 500 Gene mit der höchsten Varianz ausgewählt und erneut hierarchisch geclustert.

![Clustering der 500 variabelsten Gene](images/clustering_top500.png)

Auch bei Beschränkung auf die 500 variabelsten Gene wurden die behandelten und unbehandelten Proben vollständig voneinander getrennt.

Zur weiteren Sensitivitätsanalyse wurde das Clustering für unterschiedliche Anzahlen hochvariabler Gene wiederholt.

| Anzahl variabler Gene | Vollständige Trennung nach Behandlung |
|---:|:---:|
| 100 | nein |
| 250 | ja |
| 500 | ja |
| 1.000 | ja |
| 5.000 | ja |

Bei einer starken Reduktion auf nur 100 Gene ging die eindeutige Trennung nach Behandlungsstatus verloren. In diesem Fall wurden die behandelte und unbehandelte Probe der Zelllinie N080611 gemeinsam von den übrigen sechs Proben abgegrenzt. Ab 250 hochvariablen Genen blieb die Trennung nach Behandlungsstatus erhalten.

## Sample-Distanzen

Zur Darstellung der paarweisen Unterschiede zwischen den acht Proben wurde die euklidische Distanzmatrix zusätzlich als Heatmap visualisiert.

![Sample-Distanzen der Airway-Proben](images/sample_distance_heatmap.png)

Die Heatmap bestätigt die im Dendrogramm beobachtete Struktur. Die vier unbehandelten und die vier behandelten Proben bilden jeweils zusammenhängende Bereiche mit vergleichsweise geringen Distanzen innerhalb der jeweiligen Gruppe. Zwischen behandelten und unbehandelten Proben treten dagegen überwiegend größere Distanzen auf.

## Zusammenfassung

Das unüberwachte Clustering zeigt, dass der Behandlungsstatus mit Dexamethason eine dominante Struktur in den untersuchten Genexpressionsdaten darstellt. Bei Verwendung aller 14.224 Gene wurden die behandelten und unbehandelten Proben vollständig voneinander getrennt. Dieses Ergebnis blieb bei Complete, Average und Ward.D2 Linkage sowie bei Verwendung von mindestens 250 hochvariablen Genen erhalten.

Bei einer Reduktion auf lediglich 100 hochvariable Gene trat dagegen ein stärkerer zelllinienspezifischer Effekt hervor. Die Ergebnisse verdeutlichen damit sowohl die robuste Dexamethason-assoziierte Struktur der Genexpressionsdaten als auch den Einfluss der Feature-Auswahl auf unüberwachte Clusteranalysen.
