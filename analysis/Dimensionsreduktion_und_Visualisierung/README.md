Der Code sowie die README_Dimensionsreduktion wurden mithilfe von OpenAI o4-mini (TU BS Chatbot) erstellt.


README
======================================================================================

Hintergrund zum Datensatz
-----------------------------
Der verwendete `airway`-Datensatz basiert auf dem PLOS-ONE-Artikel
„RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells“
(doi:10.1371/journal.pone.0099625).

Der Datensatz umfasst vier primäre humane Airway-Smooth-Muscle-Zelllinien, die jeweils unbehandelt und nach Behandlung mit Dexamethason untersucht wurden. Die in diesem Skript verwendeten Verfahren PCA, t-SNE und UMAP stellen zusätzliche explorative Analysen des Datensatzes dar und reproduzieren keine Dimensionsreduktionsanalyse aus der ursprünglichen Publikation.


Zweck dieses Skripts
-----------------------
Dieses Skript liest eine VST-transformierte Gen-×-Sample-Matrix aus der CSV-Datei `VST_transformierte_Daten.csv` ein und führt darauf drei Verfahren zur Dimensionsreduktion durch:

- Lineare Projektion (PCA)
- Nicht-lineare Einbettung (t-SNE)
- Nicht-lineare Einbettung (UMAP)

Ziel ist es, die acht Proben (SRR-IDs) in 2D-Plots darzustellen und explorativ zu untersuchen, ob sich Muster und Ähnlichkeiten zwischen den Expressionsprofilen erkennen lassen.


Ablauf im Skript
-------------------

Das Analyse-Skript gliedert sich in fünf Hauptabschnitte, die im Folgenden näher erläutert werden:

1. Pakete laden  
   Zu Beginn werden die benötigten R-Bibliotheken geladen:  
   • **readr** zum Einlesen der CSV-Datei,  
   • **dplyr** für Data-Wrangling und Tabellenmanipulation,  
   • **ggplot2** für die Erstellung der Plots,  
   • **Rtsne** für t-SNE und  
   • **uwot** für UMAP.  

2. Datenimport und -strukturierung  
   Die VST-transformierten Werte werden aus `VST_transformierte_Daten.csv` eingelesen. Die ursprüngliche Tabelle enthält in der ersten Spalte `gene_id` und in den weiteren Spalten numerische Expressionswerte für die einzelnen Proben. Anschließend wird die Datenmatrix transponiert, sodass die Zeilen den Samples und die Spalten den Genen entsprechen.

3. Feature-Skalierung  
   Im vorliegenden Skript werden die Genfeatures mit `scale()` z-standardisiert. Dabei wird jede Gen-Spalte auf einen Mittelwert von 0 zentriert und auf eine Standardabweichung von 1 skaliert. Dadurch tragen Gene unabhängig von ihrer ursprünglichen Streuung auf vergleichbarer Skala zu den anschließenden Verfahren bei.

4. Dimensionsreduktion  
   • **PCA** (Principal Component Analysis) liefert eine lineare Zerlegung der Daten in Hauptkomponenten, die sukzessive möglichst viel Varianz erklären.  
   • **t-SNE** (t-Distributed Stochastic Neighbor Embedding) erzeugt eine nicht-lineare 2D-Darstellung mit Schwerpunkt auf lokalen Nachbarschaften. Die Perplexity wird im Skript an die Stichprobengröße angepasst; bei den acht vorliegenden Proben ergibt sich eine Perplexity von 2. Durch `set.seed(42)` ist die t-SNE-Darstellung reproduzierbar.  
   • **UMAP** (Uniform Manifold Approximation and Projection) erzeugt ebenfalls eine nicht-lineare 2D-Darstellung. `n_neighbors` wird mit `min(15, n_samp - 1)` an die Stichprobengröße angepasst und beträgt für die acht Proben 7. Im vorliegenden Skript wird vor UMAP kein Zufalls-Seed gesetzt, sodass die konkrete Anordnung zwischen verschiedenen Läufen variieren kann.

5. Visualisierung  
   Für jede Methode wird ein Scatterplot in zwei Dimensionen erstellt (PC1 vs. PC2, t-SNE1 vs. t-SNE2, UMAP1 vs. UMAP2). Jeder Punkt repräsentiert ein Sample und ist mit seiner SRR-ID beschriftet. Informationen zu Zelllinie und Behandlung werden in diesem Skript nicht zusätzlich über Farbe oder Symbol kodiert und müssen bei der Interpretation über die bekannten Sample-Metadaten zugeordnet werden.


Interpretation der Ergebnisse
-----------------------------

PCA
---
Im PCA-Plot werden die acht Proben auf die ersten beiden Hauptkomponenten projiziert. Bei der vorliegenden Analyse zeigt PC1 vor allem eine Trennung zwischen Dexamethason-behandelten und unbehandelten Proben, während PC2 stärker Unterschiede zwischen den Zelllinien abbildet.

Damit weist die PCA darauf hin, dass sowohl der Behandlungsstatus als auch die Zelllinienidentität zur Variation der Expressionsprofile beitragen. Die behandelten und unbehandelten Proben derselben Zelllinie sind dabei keine technischen Replikate. Ihre Abstände können daher nicht als Maß für technische Reproduzierbarkeit interpretiert werden.

Fazit zur PCA  
– PC1 zeigt vor allem eine Dexamethason-assoziierte Trennung.  
– PC2 bildet stärker zelllinienspezifische Unterschiede ab.  
– Beide Einflussfaktoren sind in den ersten beiden Hauptkomponenten sichtbar.


t-SNE
------
Der t-SNE-Plot zeigt die acht Proben in einem zweidimensionalen Raum. t-SNE betont vor allem lokale Nachbarschaften, also welche Proben im hochdimensionalen Ausgangsraum vergleichsweise ähnlich sind.

Die Achsen t-SNE1 und t-SNE2 besitzen keine direkte biologische Bedeutung. Kurze Abstände können auf lokale Ähnlichkeit hinweisen; große Abstände zwischen Punkten oder Gruppen sollten dagegen nicht als quantitatives Maß für die Stärke biologischer Unterschiede interpretiert werden.

Aufgrund der sehr kleinen Stichprobengröße von acht Proben und der entsprechend niedrigen Perplexity von 2 ist die Darstellung ausschließlich explorativ zu verstehen. Aussagen zur technischen Reproduzierbarkeit oder zur Stärke eines Treatment-Effekts lassen sich aus der t-SNE-Darstellung allein nicht ableiten.

Fazit zu t-SNE  
t-SNE ergänzt die PCA durch eine nicht-lineare Darstellung lokaler Nachbarschaften. Die konkrete Gruppierung kann Hinweise auf ähnliche Expressionsprofile liefern, sollte bei nur acht Proben jedoch vorsichtig interpretiert werden.


UMAP
----
Der UMAP-Plot projiziert die acht Proben ebenfalls in zwei Dimensionen. UMAP versucht insbesondere lokale Nachbarschaftsstrukturen abzubilden; die Achsen UMAP1 und UMAP2 besitzen dabei keine direkte biologische Bedeutung.

Nahe beieinanderliegende Punkte können auf ähnliche Expressionsprofile hinweisen. Die Abstände zwischen weiter entfernten Punkten oder Gruppen sind jedoch nicht als direktes quantitatives Maß der biologischen Unterschiede zu verstehen. Insbesondere lässt sich aus der UMAP-Darstellung nicht ableiten, dass größere Abstände automatisch größere Expressionsunterschiede bedeuten.

Bei den acht vorliegenden Proben wird `n_neighbors = 7` verwendet. Dadurch bezieht die Nachbarschaftsdefinition nahezu den gesamten Datensatz ein. Zusammen mit der sehr kleinen Stichprobengröße sollte die UMAP-Darstellung daher ausschließlich explorativ interpretiert werden.

Da im Skript vor UMAP kein `set.seed()` gesetzt wird, kann sich die konkrete Lage der Punkte bei erneutem Ausführen verändern. Aussagen über feste Richtungen entlang UMAP1 oder UMAP2 sollten deshalb vermieden werden.

Fazit zu UMAP  
UMAP liefert eine zusätzliche nicht-lineare Darstellung der Ähnlichkeiten zwischen den Proben. Aussagen über lokale Strukturen sind möglich, globale Abstände, Achsenrichtungen oder die Stärke biologischer Effekte sollten jedoch nicht überinterpretiert werden.


Gesamtfazit
-----------
PCA, t-SNE und UMAP liefern drei unterschiedliche explorative Perspektiven auf die VST-transformierten und anschließend z-standardisierten Expressionsdaten. Die PCA zeigt eine deutliche Dexamethason-assoziierte Struktur sowie zelllinienspezifische Variation. t-SNE und UMAP ergänzen diese lineare Darstellung um nicht-lineare Projektionen lokaler Nachbarschaften.

Aufgrund der sehr kleinen Stichprobengröße von acht Proben sollten insbesondere t-SNE und UMAP zurückhaltend interpretiert werden. Die Ergebnisse dienen vor allem der Visualisierung und explorativen Beschreibung der Datenstruktur.
