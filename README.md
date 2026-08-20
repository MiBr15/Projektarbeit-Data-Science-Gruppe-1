# ReadMe - Data Wrangling (AnalyseDatenstruktur.R)

## Allgemein:
Der in diesem Branch vorliegende Code "AnalyseDatenstruktur.R" dient der Untersuchung der im Bioconductor-Paket hinterlegten Rohdaten der Abhandlung und der anschließenden Umwandlung in ein für die weitere Analyse vorteilhafteres Format. Die einzelnen Schritte und Anwendungshinweise werden im Folgenden erklärt.

## Struktur und Funktion des Codes:

### Bioconductor-Paket erhalten und Paketstruktur analysieren

1. Bioconductor, Paket und Datensatz laden
2. Datenstruktur verstehen

In den genannten Codeabschnitten wurde zunächst die Struktur des Inhalts des "Airway"-Datensatzes verstanden und anschließend die Messdaten, Proben- und Geninformationen betrachtet.

### Zerlegung des Bioconductor-Pakets

3. Datenextraktion

Nach erfolgter erster Betrachtung des Datensatzes wurden die einzelnen Informationsbereiche (Messdaten, Proben-, Geninformationen) als Variablen gespeichert.

- Messwerte_roh: Messdaten des Versuchs
- Probeninformationen: Anzahl und Kennung der Zelllinien (Sample)
- Geninformationen: Informationen zu Position, Art und Namen des Gens anhand der gene_id

### Untersuchung der Struktur und Art der Verknüpfung der im Bioconductor-Format enthaltenen Daten

4. Untersuchung Messwerte_roh, Probeninformationen, Geninformationen
5. Überprüfung gleicher Sample-IDs

Zum Zusammenführen der 3 Informationsbereiche wurden in den genannten Codeabschnitten die enthaltenen Daten untersucht und ein nutzbarer Zusammenhang zwischen den 3 Informationsbereichen bestimmt. Hierbei konnte festgestellt werden, dass alle Informationsbereiche die zugehörige gene_id enthalten.

### Erstellung eines zusammenfassenden Datenformats für die weitere Analyse

6. Erstellen einer Probeninformationstabelle
7. Überprüfung erstellter Tabelle
8. Messwerte_roh ins Long Format bringen
9. Überprüfung der Tabelle
10. Messwerte_roh_long mit Tab_Probeninformationen verbinden
11. Neue Tabelle überprüfen
12. Join von Daten_roh_long und Geninformationen

Anhand der in allen Informationsbereichen vorliegenden gene_id wurden in den genannten Codeabschnitten schrittweise die 3 Informationsbereichen in eine Tabelle (long-Format) zusammengeführt. Hierfür wurden die Messwerte in eine Tabelle im Long-Format überführt und nachfolgend mit den Tabellen zu Proben- und Geninformationen anhand der gene_id passend aneinander angefügt.

### Kontrolle des neu erstellten Datensatzes für die weiterführende Analyse

13. Überprüfung des Daten_roh_long-Datensatzes

Zur Überprüfung der Richtigkeit des neu erstellten Datenformats wurde dessen struktureller Aufbau, Proben- und Genanzahl, Datenverteilung und Vorkommen von NA-Werte näher untersucht.

### Sicherung des erstellten Datenformats als csv-Datei

Nach erfolgter Prüfung des Datensatzes wurde dieser als csv-Datei gespeichert und zur weiteren Benutzung neben dem Code zur Verfügung gestellt.
