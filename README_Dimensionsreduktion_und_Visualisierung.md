Der Code sowie die README_Dimensionsreduktion wurden mithilfe von Open AI o4-mini (TU BS Chatbot) erstellt.


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

-------------------
Interpretation der Ergebnisse

Im PCA-Plot sind je vier Zelllinien in jeweils unbehandeltem (untrt) und behandeltem (trt) Zustand dargestellt. Die acht Punkte („SRR1039508“ bis „SRR1039521“) wurden auf die ersten beiden Hauptkomponenten (PC1, PC2) projiziert.

1. Hauptquelle der Varianz – PC1  
   - Entlang der horizontalen Achse (PC1) liegen die Proben der vier Zelllinien in vier deutlich separierten Clustern.  
   - Das bedeutet: Die **unterschiedlichen Zelllinien** selbst (unabhängig vom Behandlungsstatus) tragen am meisten zur Gesamtvarianz bei.  
   - PC1 kann somit als „Zelllinien‐Achse“ interpretiert werden und repräsentiert ein Ensemble von Genexpressionsmustern, die zwischen den Linien variieren (z. B. Grundexpression, angeborene Signalwege usw.).

2. Sekundäre Quelle der Varianz – PC2  
   - Die vertikale Achse (PC2) trennt innerhalb jedes Zelllinien-Clusters unbehandelte von behandelten Proben:  
     – In jedem Paar hat die trt-Probe einen höheren PC2-Wert als die untrt-Probe.  
   - Das spricht dafür, dass PC2 vor allem **die durch das Treatment induzierten Transkriptom-Veränderungen** abbildet.  
   - PC2 ist also eine „Behandlungs-Achse“: Je stärker die Proben in PC2 differieren, desto größer der Behandlungseffekt auf das Genexpressionsprofil.

3. Reproduzierbarkeit und technische Präzision  
   - Die unmittelbare Nachbarschaft je eines untrt– und trt-Paares (enger Abstand in PC1) weist auf **gute technische Replikatreue** hin.  
   - Die intra-Paar-Varianz (Abstand zwischen untrt vs. trt) ist deutlich kleiner als die inter-Zelllinien-Varianz (Abstand zwischen Clustern auf PC1).

4. Biologische Implikationen  
   - Zelllinien unterscheiden sich grundlegend im Genexpressions-„Baseline“-Profil (PC1).  
   - Die Behandlung verschiebt das Transkriptom jeder Zelllinie **in gleiche Richtung** entlang PC2, was auf einen vergleichbaren molekularen Effekt hindeutet (z. B. Aktivierung oder Repression gemeinsamer Signalwege).  
   - Trotz ähnlicher Behandlungsantwort (gleiche Richtung auf PC2) unterscheiden sich die Zelllinien weiterhin in ihrer Ausprägung (verschobene Clusterpositionen auf PC1).
Fazit  
– PC1: Hauptunterschiede zwischen den Zelllinien  
– PC2: Konsistenter Behandlungseffekt über alle Linien hinweg  
– Hohe Replikatreue innerhalb von untrt–/trt-Paaren und klare Trennung beider Effekte zeigt die Eignung der PCA zur **schnellen** Übersicht über biologische und experimentelle Einflussfaktoren.



Der t-SNE-Plot zeigt acht Proben (jeweils unbehandelt = untrt und behandelt = trt für vier Zelllinien) im zweidimensionalen t-SNE-Raum. Im Gegensatz zur PCA betont t-SNE vor allem **lokale Nachbarschaften** (also welche Proben sehr ähnlich sind), ohne die Achsen global zu interpretieren. Folgende Punkte sind bei der Auswertung relevant:

1. Bildung lokaler Cluster  
   - Die Proben jedes Zelllinien-Paares (untrt / trt) gruppieren sich eng beieinander.  
   - Das zeigt, dass innerhalb einer Zelllinie die technische Reproduzierbarkeit hoch ist und der globale Zelllinien-Unterschied gegenüber dem Behandlungs­effekt dominanter bleibt.  

2. Behandlungseffekt als kleine Verschiebung  
   - Innerhalb jedes solchen Paares liegt die behandel­te (trt) Probe leicht versetzt von der unbehandel­ten (untrt).  
   - Diese Versetzung ist in t-SNE meist nicht exakt in einer Achse abzulesen, sondern als **Richtung und Abstand** im lokalen Cluster zu interpretieren:  
     – Ein konsistentes Verschieben aller trt-Punkte weg von den untrt-Punkten weist auf einen einheitlichen transkriptionellen Response-Stil hin.  
     – Unterschiede in der Verschiebungslänge oder -richtung zwischen den Zelllinien deuten auf **zelllinien-spezifische Treatment-Antworten** hin.

3. Zelllinien-Spezifität vs. Treatment-Effekt  
   - Weil sich die vier Paare klar separiert zeigen, überlagert die **zelllinsenspezifische Baseline** die Behandlungssignale.  
   - Ein rein auf Treatment basierendes Clustering („alle trt beieinander, alle untrt beieinander“) wäre nur zu erwarten, wenn der Behandlungs­effekt die zwischen-Zelllinien-Variabilität übersteigen würde.

4. Bedeutung der Abstände  
   - **Kurze Distanzen** im t-SNE lassen auf hohe Ähnlichkeit im Original-Genexpressionsraum schließen.  
   - **Große Distanzen** zwischen Clustern (z. B. zwischen den Zelllinien-Gruppen) zeigen fundamentale Unterschiede in den Expressionsprofilen.  
   - Die exakten Achsenwerte t-SNE1/t-SNE2 haben keine unmittelbare biophysikalische Deutung, dienen nur der Visualisierung.

5. Biologische Schlussfolgerungen  
   - Die Daten bestätigen:  
     • Jede Zelllinie besitzt ein charakteristisches Expressionsprofil (vier separierte Cluster).  
     • Die Behandlung induziert in allen Zelllinien eine erkennbare, aber vergleichsweise geringere Umprogrammierung.  
     • Die **Richtung** und **Stärke** dieser Umprogrammierung können sich von Zelllinie zu Zelllinie unterscheiden (hinweisend auf unterschiedliche Sensitivitäten oder Aktivierungspfade).

Fazit  
Der t-SNE-Plot unterstreicht die starke Zelllinien-Spezifität und zeigt gleichzeitig einen konsistenten, aber zelllinien-abhängigen Behandlungs­effekt. Die enge Paarbildung (untrt/trt) belegt dabei gute Replikatreue und lokale Ähnlichkeit der Proben.




Der UMAP-Plot projiziert die acht Proben (jeweils unbehandelt = untrt und behandelt = trt für vier Zelllinien) in zwei Dimensionen, wobei sowohl lokale Nachbarschaften als auch globale Abstände im ursprünglichen hoch­dimensionalen Expressionsraum weitgehend erhalten bleiben.

 1. Globale Struktur – Zelllinienidentität  
  - Entlang der UMAP 1-Achse bilden sich vier klar separierte Gruppen, die den **jeweiligen Zelllinien-Baselines** entsprechen.  
  - Die Abstände zwischen diesen Zelllinien-Clustern sind im UMAP-Raum „echt“ vergleichbar: Je weiter zwei Cluster auseinanderliegen, desto größer sind ihre  Gesamtunterschiede im Genexpressionsprofil.  
  - Das deutet darauf hin, dass die genomische Grundverschiedenheit der Zelllinien den stärksten Beitrag zur Varianz liefert.

 2. Lokale Struktur – Treatment-Effekt  
  - Innerhalb jedes Zelllinien-Clusters liegen die unbehandelten und behandelten Proben als dichtes Paar: hohe technische Reproduzierbarkeit.  
  - Die behandelten Proben (trt) sind konsistent in eine ähnliche Richtung von den unbehandelten (untrt) verschoben – hauptsächlich entlang UMAP 2, manchmal zu einem kleineren Teil auch entlang UMAP 1.  
  - Diese Verschiebung zeigt den **transkriptomischen Effekt** der Behandlung in jeder Zelllinie an.  

 3. Zelllinien-spezifische Unterschiede in der Treatment-Antwort  
  - Länge und Richtung des Versatzes zwischen untrt und trt variieren leicht zwischen den Zelllinien:  
  - Manche Linien zeigen eine stärkere Verschiebung (größere Änderung), andere eine schwächere (moderaterer Treatment-Response).  
  - Unterschiede in der Versatz­richtung deuten an, dass in verschiedenen Zelllinien **unterschiedliche Gene** oder **Signalwege** auf das Treatment reagieren.

 4. Vergleich zu PCA und t-SNE  
  - UMAP bewahrt im Gegensatz zu t-SNE besser die **globale Geometrie**: Die relativen Abstände zwischen den Zelllinien-Clustern stimmen besser mit den Varianz­achsen der PCA überein.  
  - Gleichzeitig stellt UMAP wie t-SNE sicher, dass nahe beieinanderliegende Proben im Originalraum auch lokal nahe bleiben.

 5. Biologische Schlussfolgerungen  
  - Die **Haupttrennlinie** ist die Zelllinienidentität (UMAP 1), gefolgt von der **Behandlung** (UMAP 2).  
  - Ein konsistenter, jedoch zelllinienabhängiger Treatment-Effekt legt nahe, dass das untersuchte Präparat in allen Zelllinien wirkt, aber in Ausmaß und Mechanismus variieren kann.  
  - Keine Probe liegt isoliert fernab: Es gibt keine starken Ausreißer oder Batch-Artefakte.


---------------------------  