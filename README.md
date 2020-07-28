# bachelorarbeit
Repository für alle Programme und Dateien für die Bachelorarbeit

Infos zu den inhalten des repo:
Messungen wurden mit ising.c, skalierung.c, messfunktionen.c/h, sweeps.c/h, auswertungsfunktionen.c/h in serieller Ausführung und mit OpenMP gemacht.
Messungen mit MPI mit mpiising.c, skalierungmpi.c, mpifunktionen.c/h gemacht. 
Zusaetzliche Auswertung mit minline.c, maxlineskalierung.c.
Andere .c/h Dateien sind ältere Versionen, die zur Fehlersuche benoetigt wurde, oder Ausprobieren von Konzepten.
Programme wurden mit makefile erstellt. 
Alle Dateien mit "batch" im Titel sowie skalmpischleife.sh dienen dem Starten von Jobs auf QBiG.
Alle .gp-Dateien sind zum Plotten von Ergebnissen, dort Ausgabe in pdfs.
Im Ordner Messungen sind Messergebnisse:
	-Im Ordner Bootstrapges sind die Ergebnisse der Messungen der Observablen, Observable, Laenge und MEssungen im Titel. Dateien mit "proz" im Titel sind Messergebnisse von MPI-Messungen, alles andere von seriellen oder OpenMP-Messungen. Erzeugt in mpiising.c oder ising.c.
	-Im Ordner Mittelwerte sind die naiven Mittelwerte der Observablen gespeichert. Erzeugt in mpiising.c oder ising.c.
	-Im Ordner Zeiten sind Messungen, wie lange die Programmausführung verschiedener Programm(teile) gebraucht hat. Die meisten wurden mit skalierung.c oder skalierungmpi.c erstellt, einige auch von ising.c.
Im Ordner outanderror sind .out-Dateien enthalten, die bei Jobs auf QBiG erzeugt wurden.
Im Ordner texdateien sind die tex-Dokumente, aus denen die Bachelorarbeit erzeugt wird, sowie einige zusätzliche Hilfsdokumente, wie die Datei mit den Literatureinträgen oder die Pakete von Ian Brock, mit denen die Arbeit gesetzt wurde.
	-Im Ordner Bilder sind die Dateien plot.gp und multiplot.gp, mit denen die Bilder für die Arbeit erzeugt wurden, sowie die Messwerte in .txt-Dateien, die dafuer benoetigt wurden.

Infos zum Ausführen:	
Zum Ausführen der Programme müssen im Ordner Messungen die Unterordner MPIMessungen, ThermalisierteGitter und Messwerte erstellt werden und die Programme, z.B. mit dem makefile, kompiliert werden.
Folgende Parameter sollten bei der Ausführung übergeben werden, in dieser reihenfolge:
ising: Gitterlänge, Anzahl an verfügbaren OpenMP-Threads, Schrittweite der Temperaturen
skalierung: Gitterlänge, Temperatur, ein Merkmal, das in den Dateinamen der Messungen geschrieben wird
mpiising: Gitterlänge, Schrittweite der Temperaturen
mpiskalierung: Gitterlänge, Temperatur, ein Merkmal, das in den Dateinamen der Messungen geschrieben wird
Zum Erstellen der Bachelorarbeit müssen die Programme plot.gp und multiplot.gp ausgeführt werden und die Datei Bachelorarbeit.tex kompiliert werden.
