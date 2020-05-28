set ter pdfcairo size 5 in, 3.5 in

set out "Messungen/speedup.pdf"
set xlabel "Anzahl Prozessoren"
set ylabel "Zeit(n Cores)/Zeit(1 Core)"

set title "Laenge=50"
plot "Messungen/Zeiten/zeitenmessen-laenge0050-m100000.txt" u 1:4 lt 7 title "" 
set title "Laenge=100"
plot "Messungen/Zeiten/zeitenmessen-laenge0100-m010000-mehrere.txt" u 1:4 lt 7 title "" 

set ylabel "Zeit(1 Core)/Zeit(n Cores)"

set title "Laenge=50"
plot "Messungen/Zeiten/zeitenmessen-laenge0050-m100000.txt" u 1:(1/$4) lt 7 title "" 
set title "Laenge=100"
plot "Messungen/Zeiten/zeitenmessen-laenge0100-m010000-mehrere.txt" u 1:(1/$4) lt 7 title "" 


set out "Messungen/Zeiten/zeitprocore.pdf"

set style line 1 lc 7 lt 7 
set style line 2 lc 1 lt 7
set style line 3 lc 5 lt 7


set xlabel "Anzahl cores"
set ylabel "Benötigte Zeit"

plot "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-01.txt" u (($1==4)?1:1/0):4 ls 1 title "Gesamtzeit",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-02.txt" u (($1==4)?2:1/0):4 ls 1 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-03.txt" u (($1==4)?3:1/0):4 ls 1 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-04.txt" u (($1==4)?4:1/0):4 ls 1 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-05.txt" u (($1==4)?5:1/0):4 ls 1 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-06.txt" u (($1==4)?6:1/0):4 ls 1 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-07.txt" u (($1==4)?7:1/0):4 ls 1 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-08.txt" u (($1==4)?8:1/0):4 ls 1 title "",\
     "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-01.txt" u (($1==3)?1:1/0):4 ls 2 title "Messen",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-02.txt" u (($1==3)?2:1/0):4 ls 2 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-03.txt" u (($1==3)?3:1/0):4 ls 2 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-04.txt" u (($1==3)?4:1/0):4 ls 2 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-05.txt" u (($1==3)?5:1/0):4 ls 2 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-06.txt" u (($1==3)?6:1/0):4 ls 2 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-07.txt" u (($1==3)?7:1/0):4 ls 2 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-08.txt" u (($1==3)?8:1/0):4 ls 2 title "",\
     "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-01.txt" u (($1==2)?1:1/0):4 ls 3 title "Bootstrappen",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-02.txt" u (($1==2)?2:1/0):4 ls 3 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-03.txt" u (($1==2)?3:1/0):4 ls 3 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-04.txt" u (($1==2)?4:1/0):4 ls 3 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-05.txt" u (($1==2)?5:1/0):4 ls 3 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-06.txt" u (($1==2)?6:1/0):4 ls 3 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-07.txt" u (($1==2)?7:1/0):4 ls 3 title "",\
	 "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-08.txt" u (($1==2)?8:1/0):4 ls 3 title "",\
	 
set title "1 core"
set xlabel "Temperatur"
plot "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-01.txt" u (($1==0)?$2:1/0):4 ls 1 title "Messen 10.000"
plot "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-01.txt" u (($1==1)?$2:1/0):4 ls 1 title "Bootstrappen 10.000"
set xlabel "Temperatur"
plot "Messungen/Zeiten/zeiten-laenge-0050-m-000100-cores-01.txt" u (($1==0)?$2:1/0):4 ls 1 title "Messen 100"
plot "Messungen/Zeiten/zeiten-laenge-0050-m-000100-cores-01.txt" u (($1==1)?$2:1/0):4 ls 1 title "Bootstrappen 100"
set xlabel "Temperatur"
plot "Messungen/Zeiten/zeiten-laenge-0050-m-001000-cores-01.txt" u (($1==0)?$2:1/0):4 ls 1 title "Messen 1000"
plot "Messungen/Zeiten/zeiten-laenge-0050-m-001000-cores-01.txt" u (($1==1)?$2:1/0):4 ls 1 title "Bootstrappen 1000"

set title "2 cores"
set xlabel "Temperatur"
plot "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-02.txt" u (($1==0)?$2:1/0):4 ls 1 title "Messen"
plot "Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-02.txt" u (($1==1)?$2:1/0):4 ls 1 title "Bootstrappen"

#plot "Messungen/Zeiten/zeitenalle-laenge-0050-m-010000.txt" u 1:6 title "Gesamtzeit"

set out "Messungen/Zeiten/zeitblocksgenerieren.pdf"
do for [zahl=0:20:5]{
datei = sprintf("Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-02-%.2d.txt", zahl)
titel=sprintf("Ohne Messungen vorher %d", zahl/5)
set title titel
plot datei u (($1==1)?$2:1/0):4 ls 1 title ""
}
do for [zahl=1:21:5]{
datei = sprintf("Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-02-%.2d.txt", zahl)
titel=sprintf("Mit Mesungen vorher %d", (zahl-1)/5)
set title titel
plot datei u (($1==1)?$2:1/0):4 ls 1 title ""
}
do for [zahl=2:17:5]{
datei = sprintf("Messungen/Zeiten/zeiten-laenge-0050-m-010000-cores-02-%.2d.txt", zahl)
titel=sprintf("Mit Mesungen vorher %d", (zahl-2)/5)
set title titel
plot datei u (($1==1)?$2:1/0):4 ls 1 title ""
}
