# Clustering-Analyse

## Fragestellung

Mittels unüberwachtem Clustering wurde untersucht, ob sich die acht Airway-Proben allein anhand ihrer Genexpressionsprofile in biologisch sinnvolle Gruppen einteilen lassen. Im Mittelpunkt stand die Frage, ob die Clusterstruktur eher dem Dexamethason-Behandlungsstatus oder den vier untersuchten Zelllinien entspricht.

## Hypothese

Die Behandlung mit Dexamethason führt zu systematischen Veränderungen des Genexpressionsprofils, sodass behandelte und unbehandelte Proben auch ohne Verwendung der bekannten Gruppeninformation getrennt werden können.

## Verwendete Daten

Für die Analyse werden folgende Dateien benötigt:

- `VST_transformierte_Daten.csv`
- `Probeninformationen_Normalisierung.csv`

Die Expressionsmatrix enthält 14.224 nach der Qualitätskontrolle verbliebene Gene für acht Proben.

## Vorgehen

1. Berechnung euklidischer Distanzen zwischen den acht Proben auf Basis der VST-transformierten Expressionswerte.
2. Hierarchisches Clustering mit Complete Linkage.
3. Robustheitsprüfung mit Average Linkage und Ward.D2.
4. Auswahl hochvariabler Gene anhand der Varianz über alle acht Proben.
5. Wiederholung des Clusterings mit den 100, 250, 500, 1.000 und 5.000 variabelsten Genen.
6. Visualisierung der vollständigen Sample-Distanzmatrix als Heatmap.

## Ergebnisse

Bei Verwendung aller 14.224 Gene trennt das hierarchische Clustering die Proben vollständig nach Behandlungsstatus. Alle vier unbehandelten Proben bilden einen Hauptcluster, während alle vier mit Dexamethason behandelten Proben dem zweiten Hauptcluster zugeordnet werden.

Die vollständige Trennung bleibt bei Complete, Average und Ward.D2 Linkage erhalten und ist damit gegenüber der Wahl des Linkage-Verfahrens robust.

Auch bei Verwendung der 250, 500, 1.000 und 5.000 variabelsten Gene ergibt sich eine vollständige Trennung nach Behandlungsstatus. Bei einer starken Reduktion auf die 100 variabelsten Gene geht diese eindeutige Treatment-Trennung dagegen verloren. In diesem Fall werden die behandelte und unbehandelte Probe der Zelllinie N080611 gemeinsam von den übrigen sechs Proben abgegrenzt.

Die Analyse zeigt damit eine robuste Dexamethason-assoziierte Struktur der Genexpressionsdaten, verdeutlicht jedoch gleichzeitig den Einfluss einer sehr starken Feature-Reduktion auf das Ergebnis eines unüberwachten Clusterings.

## R-Pakete

- `readr`
- `dplyr`
- `tibble`
- `pheatmap`
