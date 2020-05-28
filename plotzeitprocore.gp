set ter pdfcairo size 5 in, 3.5 in

set xlabel "Anzahl cores"
set ylabel "Benötigte Zeit"

set out "Messungen/Zeitprocore.pdf"

set style line 1 lc 7 lt 7 
set style line 2 lc 1 lt 7
set style line 3 lc 5 lt 7

plot "Messungen/zeiten-laenge-0050-m-010000-cores-01.txt" u (($1==4)?1:1/0):4 ls 1 title "Gesamtzeit",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-02.txt" u (($1==4)?2:1/0):4 ls 1 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-03.txt" u (($1==4)?3:1/0):4 ls 1 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-04.txt" u (($1==4)?4:1/0):4 ls 1 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-05.txt" u (($1==4)?5:1/0):4 ls 1 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-06.txt" u (($1==4)?6:1/0):4 ls 1 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-07.txt" u (($1==4)?7:1/0):4 ls 1 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-08.txt" u (($1==4)?8:1/0):4 ls 1 title "",\
     "Messungen/zeiten-laenge-0050-m-010000-cores-01.txt" u (($1==3)?1:1/0):4 ls 2 title "Messen",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-02.txt" u (($1==3)?2:1/0):4 ls 2 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-03.txt" u (($1==3)?3:1/0):4 ls 2 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-04.txt" u (($1==3)?4:1/0):4 ls 2 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-05.txt" u (($1==3)?5:1/0):4 ls 2 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-06.txt" u (($1==3)?6:1/0):4 ls 2 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-07.txt" u (($1==3)?7:1/0):4 ls 2 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-08.txt" u (($1==3)?8:1/0):4 ls 2 title "",\
     "Messungen/zeiten-laenge-0050-m-010000-cores-01.txt" u (($1==2)?1:1/0):4 ls 3 title "Bootstrappen",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-02.txt" u (($1==2)?2:1/0):4 ls 3 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-03.txt" u (($1==2)?3:1/0):4 ls 3 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-04.txt" u (($1==2)?4:1/0):4 ls 3 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-05.txt" u (($1==2)?5:1/0):4 ls 3 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-06.txt" u (($1==2)?6:1/0):4 ls 3 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-07.txt" u (($1==2)?7:1/0):4 ls 3 title "",\
	 "Messungen/zeiten-laenge-0050-m-010000-cores-08.txt" u (($1==2)?8:1/0):4 ls 3 title "",\
	 
set title "1 core"
set xlabel "Temperatur"
plot "Messungen/zeiten-laenge-0050-m-010000-cores-01.txt" u (($1==0)?$2:1/0):4 ls 1 title "Messen 10.000"
plot "Messungen/zeiten-laenge-0050-m-010000-cores-01.txt" u (($1==1)?$2:1/0):4 ls 1 title "Bootstrappen 10.000"
set xlabel "Temperatur"
plot "Messungen/zeiten-laenge-0050-m-000100-cores-01.txt" u (($1==0)?$2:1/0):4 ls 1 title "Messen 100"
plot "Messungen/zeiten-laenge-0050-m-000100-cores-01.txt" u (($1==1)?$2:1/0):4 ls 1 title "Bootstrappen 100"
set xlabel "Temperatur"
plot "Messungen/zeiten-laenge-0050-m-001000-cores-01.txt" u (($1==0)?$2:1/0):4 ls 1 title "Messen 1000"
plot "Messungen/zeiten-laenge-0050-m-001000-cores-01.txt" u (($1==1)?$2:1/0):4 ls 1 title "Bootstrappen 1000"

set title "2 cores"
set xlabel "Temperatur"
plot "Messungen/zeiten-laenge-0050-m-010000-cores-02.txt" u (($1==0)?$2:1/0):4 ls 1 title "Messen"
plot "Messungen/zeiten-laenge-0050-m-010000-cores-02.txt" u (($1==1)?$2:1/0):4 ls 1 title "Bootstrappen"

#plot "Messungen/zeitenalle-laenge-0050-m-010000.txt" u 1:6 title "Gesamtzeit"

set out "Messungen/zeitblocksgenerieren.pdf"
do for [zahl=0:20:5]{
datei = sprintf("Messungen/zeiten-laenge-0050-m-010000-cores-02-%.2d.txt", zahl)
titel=sprintf("Ohne Messungen vorher %d", zahl/5)
set title titel
plot datei u (($1==1)?$2:1/0):4 ls 1 title ""
}
do for [zahl=1:21:5]{
datei = sprintf("Messungen/zeiten-laenge-0050-m-010000-cores-02-%.2d.txt", zahl)
titel=sprintf("Mit Mesungen vorher %d", (zahl-1)/5)
set title titel
plot datei u (($1==1)?$2:1/0):4 ls 1 title ""
}
do for [zahl=2:22:5]{
datei = sprintf("Messungen/zeiten-laenge-0050-m-010000-cores-02-%.2d.txt", zahl)
titel=sprintf("Mit Mesungen vorher %d", (zahl-2)/5)
set title titel
plot datei u (($1==1)?$2:1/0):4 ls 1 title ""
}
