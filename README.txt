README  
======================================================================================

Hintergrund aus dem Paper  
-----------------------------  
Im PLOS-ONE-Artikel “Single-cell RNA-seq reveals dynamic, random monoallelic gene expression in mammalian cells” (doi:10.1371/journal.pone.0099625) wenden die Autoren u. a. PCA und nicht-lineare Einbettungsverfahren an, um hochdimensionale Genexpressionsdaten zu visualisieren und Zell-Subpopulationen zu identifizieren. Sie zeigen, wie VST-transformierte RNA-Seq Counts durch Dimensionsreduktion in zwei oder drei Komponenten projiziert werden können, um Muster, Cluster und Übergangszustände im Datensatz aufzudecken.

Zweck dieses Skripts  
-----------------------  
Dieses Skript nimmt eine VST-transformierte Gen-×-Sample-Matrix (CSV „VST_transformierte_Daten.csv“) und reproduziert die wesentlichen Schritte aus dem Paper:

- Lineare Projektion (PCA)  
- Nicht-lineare Einbettung (t-SNE, UMAP)  

Ziel ist es, die Proben (z. B. SRR IDs) in 2D-Plots darzustellen, um Gruppen­bildung und Ähnlichkeiten analog zur Publikation zu explorieren.

Ablauf im Skript  
-------------------

Der Analyse-Skript gliedert sich in fünf Hauptabschnitte, die im Folgenden ausführlicher erläutert werden:

1. Pakete laden  
   Zu Beginn werden alle benötigten R-Bibliotheken geladen:  
   • **readr** zum schnellen Einlesen der CSV-Datei,  
   • **dplyr** für Data-Wrangling und Tabellenmanipulation,  
   • **ggplot2** für die Erstellung publication-quality Plots,  
   • **Rtsne** für t-SNE und  
   • **uwot** für UMAP.  
   Damit stehen effiziente Werkzeuge zur Datenvorbereitung und verschiedenen Einbettungsmethoden bereit.

2. Datenimport und -strukturierung  
   Die VST-transformierten Werte werden aus `VST_transformierte_Daten.csv` eingelesen. Die ursprüngliche Tabelle enthält in der ersten Spalte `gene_id` und in den weiteren Spalten numerische Expressionswerte für jede Probe (z. B. SRR1039508). Anschließend wird die Datenmatrix so aufgebaut, dass Zeilen den Samples entsprechen und Spalten die Gene repräsentieren (Transposition der Einlesetabelle).

3. Feature-Skalierung  
   Bevor eine Dimensionsreduktion sinnvoll durchgeführt werden kann, müssen die Gen-Expressionsprofile vergleichbar gemacht werden. Dazu wird jede Spalte (jedes Gen) mittels Z-Score-Normierung zentriert (Mittel = 0) und auf Ein-Standardabweichung skaliert. So wird verhindert, dass Gene mit hohen absoluten VST-Werten die Projektion dominieren.

4. Dimensionsreduktion  
   • **PCA** (Principal Component Analysis) liefert eine lineare Zerlegung der Daten in die beiden Hauptkomponenten, die die größte Varianz erklären.  
   • **t-SNE** (t-Distributed Stochastic Neighbor Embedding) erzeugt eine nicht-lineare 2D-Darstellung, bei der lokale Nachbarschaften erhalten bleiben. Die Perplexity wird automatisch an die Stichprobengröße angepasst (max. ≈ (n–1)/3), um Über- oder Unter-Fitting zu vermeiden.  
   • **UMAP** (Uniform Manifold Approximation and Projection) kombiniert globale und lokale Strukturinformation; auch hier wird `n_neighbors` dynamisch so gewählt, dass es kleiner als (Anzahl Samples − 1) ist.

5. Visualisierung  
   Für jede Methode wird ein Scatterplot in 2D erstellt (PC1 vs. PC2, t-SNE1 vs. t-SNE2, UMAP1 vs. UMAP2). Jeder Punkt repräsentiert ein Sample und ist mit seiner ID beschriftet. Die Anordnung der Punkte erlaubt die Identifikation von Clustern, Subpopulationen oder Ausreißern analog zu den Analysen in doi:10.1371/journal.pone.0099625.

Interpretation der Ergebnisse  
--------------------------------  

- Punkte-Cluster in PCA bestätigen, welche Proben ähnliche Varianzachsen teilen (globaler Überblick).  
- t-SNE zeigt feinere Gruppierung, oft Zell-Subtypen oder technische Chargen.  
- UMAP verbindet die Vorteile beider Welten: Erhält globale Strukturen und löst lokale Cluster gut auf.  

Analog zur PLOS-ONE-Publikation hilft diese Pipeline, aus hunderten oder tausenden Genfeatures sinnvolle 2D-Darstellungen zu erzeugen, Subpopulationen zu entdecken und Hypothesen über zelluläre Zustände abzuleiten.
