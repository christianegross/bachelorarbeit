set ter pdfcairo size 5 in, 3.5 in

set xlabel "Anzahl Prozessoren"
set ylabel "Zeit(1 Core)/Zeit(n Cores)"
set key top left

set out "Messungen/skalierungkeinbackslash.pdf"

set title "T=0,5"

plot "Messungen/Zeiten/zmittel-m001000-node02v5t05keinescpar-l12.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 6 ps 0.4 title "L=12",\
	 "Messungen/Zeiten/zmittel-m001000-node02v5t05keinescpar-l36.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 12 ps 0.4 title "L=36",\
	 "Messungen/Zeiten/zmittel-m001000-node02v5t05keinescpar-l48.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 7 ps 0.4 title "L=48",\
	 "Messungen/Zeiten/zmittel-m001000-node02v5t05keinescpar-l96.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 8 ps 0.4 title "L=96",\
	 "Messungen/Zeiten/zmittel-m001000-node02v5t05keinescpar-l768.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 11 ps 0.4 title "L=768",\
#	 "Messungen/Zeiten/zmittel-m001000-node02v5t05keinescpar-l192.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 9 ps 0.4 title "L=192",\
#	 "Messungen/Zeiten/zmittel-m001000-node02v5t05keinescpar-l384.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 10 ps 0.4 title "L=384",\

set out "Messungen/versionensweepmehreregeneratoren"

set title "L=100, T=0,5"

plot "Messungen/Zeiten/zmittel-m001000-node02sweepmgv1-l100.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 6 ps 0.4 title "v1",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv2-l100.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 7 ps 0.4 title "v2",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv3-l100.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 8 ps 0.4 title "v3",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv4-l100.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 9 ps 0.4 title "v4",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv5-l100.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 10 ps 0.4 title "v5",\


set title "L=500, T=0,5"

plot "Messungen/Zeiten/zmittel-m001000-node02sweepmgv1-l500.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 6 ps 0.4 title "v1",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv2-l500.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 7 ps 0.4 title "v2",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv3-l500.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 8 ps 0.4 title "v3",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv4-l500.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 9 ps 0.4 title "v4",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv5-l500.txt" u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 10 ps 0.4 title "v5",\

set key top right
set title "L=100, T=0,5"

set ylabel "Zeit(1000 Messungen)/s"


plot "Messungen/Zeiten/zmittel-m001000-node02sweepmgv1-l100.txt" u (($1!=0)?$1:1/0):3:4 w yerrorbars lt 6 ps 0.4 title "v1",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv2-l100.txt" u (($1!=0)?$1:1/0):3:4 w yerrorbars lt 7 ps 0.4 title "v2",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv3-l100.txt" u (($1!=0)?$1:1/0):3:4 w yerrorbars lt 8 ps 0.4 title "v3",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv4-l100.txt" u (($1!=0)?$1:1/0):3:4 w yerrorbars lt 9 ps 0.4 title "v4",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv5-l100.txt" u (($1!=0)?$1:1/0):3:4 w yerrorbars lt 10 ps 0.4 title "v5",\

set title "L=500, T=0,5"

plot "Messungen/Zeiten/zmittel-m001000-node02sweepmgv1-l500.txt" u (($1!=0)?$1:1/0):3:4 w yerrorbars lt 6 ps 0.4 title "v1",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv2-l500.txt" u (($1!=0)?$1:1/0):3:4 w yerrorbars lt 7 ps 0.4 title "v2",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv3-l500.txt" u (($1!=0)?$1:1/0):3:4 w yerrorbars lt 8 ps 0.4 title "v3",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv4-l500.txt" u (($1!=0)?$1:1/0):3:4 w yerrorbars lt 9 ps 0.4 title "v4",\
	 "Messungen/Zeiten/zmittel-m001000-node02sweepmgv5-l500.txt" u (($1!=0)?$1:1/0):3:4 w yerrorbars lt 10 ps 0.4 title "v5",\

set out "Messungen/skalierungtemperaturen.pdf"

set xlabel "Anzahl Prozessoren"
set ylabel "Zeit(1 Core)/Zeit(n Cores)"
set key top left
set title "L=500, zwei Schleifen, ohne leere Durchläufe"
plot "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t1.txt" u ((($2==500)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 6 ps 0.4 title "T=0,5",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t2.txt" u ((($2==500)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 7 ps 0.4 title "T=2,0",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t3.txt" u ((($2==500)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 8 ps 0.4 title "T=2,5",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t4.txt" u ((($2==500)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 9 ps 0.4 title "T=4,5",\

set key top right
set title "L=10, zwei Schleifen, ohne leere Durchläufe"
plot "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t1.txt" u ((($2==10)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 6 ps 0.4 title "T=0,5",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t2.txt" u ((($2==10)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 7 ps 0.4 title "T=2,0",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t3.txt" u ((($2==10)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 8 ps 0.4 title "T=2,5",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t4.txt" u ((($2==10)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 9 ps 0.4 title "T=4,5",\
	 
set title "L=500, zwei Schleifen, ohne leere Durchläufe"
set yrange [0:22]
set ylabel "Zeit(1000 Messungen)/s"

plot "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t1.txt" u ((($2==500)&&($1!=0)&&($7==0.5))?$1:1/0):3:4 w yerrorbars lt 6 ps 0.4 title "T=0,5",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t2.txt" u ((($2==500)&&($1!=0)&&($7==2.0))?$1:1/0):3:4 w yerrorbars lt 7 ps 0.4 title "T=2,0",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t3.txt" u ((($2==500)&&($1!=0)&&($7==2.5))?$1:1/0):3:4 w yerrorbars lt 8 ps 0.4 title "T=2,5",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t4.txt" u ((($2==500)&&($1!=0)&&($7==4.5))?$1:1/0):3:4 w yerrorbars lt 9 ps 0.4 title "T=4,5"
unset yrange
	 
set title "L=500 - mehrere Messmethoden"
unset yrange
set ylabel "Zeit(1 Core)/Zeit(n Cores)"
plot "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t4.txt" u ((($2==500)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 6 ps 0.4 title "zwei Schleifen ohne leere durchläufe",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-eineschleife-t4.txt" u ((($2==500)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 8 ps 0.4 title "eine Schleife",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-leeredurchlaeufe-t4.txt" u ((($2==500)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 7 ps 0.4 title "zwei Schleifen mit leeren Durchläufen"

set yrange [0:22]
set ylabel "Zeit(1000 Messungen)/s"
set title "L=500 - absolute Zeit"
plot "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t4.txt" u ((($2==500)&&($1!=0))?$1:1/0):3:4 w yerrorbars lt 6 ps 0.4 title "zwei Schleifen ohne leere durchläufe",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-eineschleife-t4.txt" u ((($2==500)&&($1!=0))?$1:1/0):3:4 w yerrorbars lt 8 ps 0.4 title "eine Schleife",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg-leeredurchlaeufe-t4.txt" u ((($2==500)&&($1!=0))?$1:1/0):3:4 w yerrorbars lt 7 ps 0.4 title "zwei Schleifen mit leeren Durchläufen"
unset yrange

set out "Messungen/skalierungtryflip.pdf"

set title "L=500, nur ziehen, wenn noetig"
plot "Messungen/Zeiten/zmittel-m001000-node02mg2snichtimmerzufall-mehreret.txt" u ((($2==500)&&($1!=0)&&($7==0.5))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 6 ps 0.4 title "T=0,5",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg2snichtimmerzufall-mehreret.txt" u ((($2==500)&&($1!=0)&&($7==2.0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 7 ps 0.4 title "T=2,0",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg2snichtimmerzufall-mehreret.txt" u ((($2==500)&&($1!=0)&&($7==2.5))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 8 ps 0.4 title "T=2,5",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg2snichtimmerzufall-mehreret.txt" u ((($2==500)&&($1!=0)&&($7==4.5))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 9 ps 0.4 title "T=4,5",\
	  "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t1.txt" u ((($2==500)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 10 ps 0.4 title "T=0,5, immer ziehen",\
	  "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t4.txt" u ((($2==500)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 11 ps 0.4 title "T=4,5, immer ziehen",\
	 
set title "L=500, nur ziehen, wenn noetig"
set key top right
plot "Messungen/Zeiten/zmittel-m001000-node02mg2snichtimmerzufall-mehreret.txt" u ((($2==500)&&($1!=0)&&($7==0.5))?$1:1/0):3:4 w yerrorbars lt 6 ps 0.4 title "T=0,5",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg2snichtimmerzufall-mehreret.txt" u ((($2==500)&&($1!=0)&&($7==2.0))?$1:1/0):3:4 w yerrorbars lt 7 ps 0.4 title "T=2,0",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg2snichtimmerzufall-mehreret.txt" u ((($2==500)&&($1!=0)&&($7==2.5))?$1:1/0):3:4 w yerrorbars lt 8 ps 0.4 title "T=2,5",\
	 "Messungen/Zeiten/zmittel-m001000-node02mg2snichtimmerzufall-mehreret.txt" u ((($2==500)&&($1!=0)&&($7==4.5))?$1:1/0):3:4 w yerrorbars lt 9 ps 0.4 title "T=4,5",\
	  "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t1.txt" u ((($2==500)&&($1!=0))?$1:1/0):3:4 w yerrorbars lt 10 ps 0.4 title "T=0,5, immer ziehen",\
	  "Messungen/Zeiten/zmittel-m001000-node02mg-inline-t4.txt" u ((($2==500)&&($1!=0))?$1:1/0):3:4 w yerrorbars lt 11 ps 0.4 title "T=4,5, immer ziehen",\

set out "Messungen/Skalierung_mehrere_generatoren.pdf"


set xlabel "Anzahl Prozessoren"
set ylabel "Zeit(n Cores)/Zeit(1 Core)"
set key top right
#set title "lcpunode01"
set title "Laenge=500"
plot "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode01mehreregeneratoren.txt" u ((($2==500)&&($1!=0))?$1:1/0):5:6 w yerrorbars lt 6 ps 0.4 title "lcpunode01",\
	 "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode02mehreregeneratoren.txt" u ((($2==500)&&($1!=0))?$1:1/0):5:6 w yerrorbars lt 7 ps 0.4 title "lcpunode02"


set ylabel "Zeit(1 Core)/Zeit(n Cores)"
set key top left
plot "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode01mehreregeneratoren.txt" u ((($2==500)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 6 ps 0.4 title "lcpunode01",\
	 "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode02mehreregeneratoren.txt" u ((($2==500)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 7 ps 0.4 title "lcpunode02"

set xlabel "Anzahl Prozessoren"
set ylabel "Zeit(n Cores)/Zeit(1 Core)"
set key top right
#set title "lcpunode01"
set title "Laenge=100"
plot "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode01mehreregeneratoren.txt" u ((($2==100)&&($1!=0))?$1:1/0):5:6 w yerrorbars lt 6 ps 0.4 title "lcpunode01",\
	 "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode02mehreregeneratoren.txt" u ((($2==100)&&($1!=0))?$1:1/0):5:6 w yerrorbars lt 7 ps 0.4 title "lcpunode02"


set ylabel "Zeit(1 Core)/Zeit(n Cores)"
set key top left
plot "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode01mehreregeneratoren.txt" u ((($2==100)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 6 ps 0.4 title "lcpunode01",\
	 "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode02mehreregeneratoren.txt" u ((($2==100)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 7 ps 0.4 title "lcpunode02"

set xlabel "Anzahl Prozessoren"
set ylabel "Zeit(n Cores)/Zeit(1 Core)"
set key top right
#set title "lcpunode01"
set title "Laenge=10"
plot "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode01mehreregeneratoren.txt" u ((($2==10)&&($1!=0))?$1:1/0):5:6 w yerrorbars lt 6 ps 0.4 title "lcpunode01",\
	 "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode02mehreregeneratoren.txt" u ((($2==10)&&($1!=0))?$1:1/0):5:6 w yerrorbars lt 7 ps 0.4 title "lcpunode02"


set ylabel "Zeit(1 Core)/Zeit(n Cores)"
set key top left
plot "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode01mehreregeneratoren.txt" u ((($2==10)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 6 ps 0.4 title "lcpunode01",\
	 "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode02mehreregeneratoren.txt" u ((($2==10)&&($1!=0))?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars lt 7 ps 0.4 title "lcpunode02"

set out "Messungen/skalierunglaenge.pdf"

f1(x)=c1*(x-a1)**2+b1
fit f1(x) "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufevmdynamic.txt" u (($1==1)?$2:1/0):3:4 yerrors via a1, b1, c1

f2(x)=c2*(x-a2)**2+b2
fit f2(x) "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode00dynamicchunk2.txt" u (($1==1)?$2:1/0):3:4 yerrors via a2, b2, c2
f3(x)=c3*(x-a3)**2+b3
fit f3(x) "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode00dynamicchunk2.txt" u (($1==2)?$2:1/0):3:4 yerrors via a3, b3, c3

set xlabel "Laenge"
set ylabel "Zeit"

plot "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufevmdynamic.txt" u (($1==1)?$2:1/0):3:4 w yerrorbars lt 7 ps 0.3 title "1 Kern", f1(x) lt 6 title "Anpassung"
plot "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode00dynamicchunk2.txt" u (($1==1)?$2:1/0):3:4 w yerrorbars lt 7 ps 0.3 title "1 Kern", f2(x) lt 6 title "Anpassung",\
	 "Messungen/Zeiten/zeitenmittel-m001000-mehrerelaengenunddurchlaeufenode00dynamicchunk2.txt" u (($1==2)?$2:1/0):3:4 w yerrorbars lt 8 ps 0.3 title "2 Kerne", f3(x) lt 6 title ""

set out "Messungen/speedupmitlaengen.pdf"
set xlabel "Anzahl Prozessoren"
set ylabel "Zeit(n Cores)/Zeit(1 Core)"

set title "Vergleich mehrere Längen"
set key top left
set xrange[0.5:9.5]
plot "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==6)?$1:1/0):4 lt 7 title "laenge=6",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==10)?$1:1/0):4 lt 1 title "laenge=10",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==20)?$1:1/0):4 lt 8 title "laenge=20",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==36)?$1:1/0):4 lt 9 title "laenge=35",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==50)?$1:1/0):4 lt 2 title "laenge=50",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==100)?$1:1/0):4 lt 3 title "laenge=100",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==200)?$1:1/0):4 lt 4 title "laenge=200",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==350)?$1:1/0):4 lt 5 title "laenge=350",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==500)?$1:1/0):4 lt 6 title "laenge=500"

set title "Vergleich mehrere Längen(Ausschnitt)"	 
set yrange[0:4]
replot 
unset yrange
unset xrange

set ylabel "Zeit(1 Cores)/Zeit(n Cores)"

set title "Vergleich mehrere Längen"
set key top left
set xrange[0:9.5]
plot "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==6)?$1:1/0):(1/$4) lt 7 title "laenge=6",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==10)?$1:1/0):(1/$4) lt 1 title "laenge=10",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==20)?$1:1/0):(1/$4) lt 8 title "laenge=20",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==36)?$1:1/0):(1/$4) lt 9 title "laenge=35",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==50)?$1:1/0):(1/$4) lt 2 title "laenge=50",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==100)?$1:1/0):(1/$4) lt 3 title "laenge=100",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==200)?$1:1/0):(1/$4) lt 4 title "laenge=200",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==350)?$1:1/0):(1/$4) lt 5 title "laenge=350",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u (($5==500)?$1:1/0):(1/$4) lt 6 title "laenge=500"

set title "lcpunode02"	 
plot "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufelcpunode02.txt" u (($5==6)?$1:1/0):(1/$4) lt 7 title "laenge=6",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufelcpunode02.txt" u (($5==10)?$1:1/0):(1/$4) lt 1 title "laenge=10",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufelcpunode02.txt" u (($5==20)?$1:1/0):(1/$4) lt 8 title "laenge=20",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufelcpunode02.txt" u (($5==36)?$1:1/0):(1/$4) lt 9 title "laenge=35",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufelcpunode02.txt" u (($5==50)?$1:1/0):(1/$4) lt 2 title "laenge=50",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufelcpunode02.txt" u (($5==100)?$1:1/0):(1/$4) lt 3 title "laenge=100",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufelcpunode02.txt" u (($5==200)?$1:1/0):(1/$4) lt 4 title "laenge=200",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufelcpunode02.txt" u (($5==350)?$1:1/0):(1/$4) lt 5 title "laenge=350",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufelcpunode02.txt" u (($5==500)?$1:1/0):(1/$4) lt 6 title "laenge=500"
	 
plot "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufelcpunode02.txt" u (($5==200)?$1:1/0):(1/$4) lt 7 title "laenge=200",\
	 "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufelcpunode02.txt" u (($5==500)?$1:1/0):(1/$4) lt 1 title "laenge=500"



set out "Messungen/speedupmitlaengengemittelt.pdf"
set title ""
plot "Messungen/Zeiten/zeitenmehrerelaengenmittelnaiv.txt" u (($2==6)?$1:1/0):(1/$3) lt 7 title "laenge=6",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaiv.txt" u (($2==10)?$1:1/0):(1/$3) lt 1 title "laenge=10",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaiv.txt" u (($2==20)?$1:1/0):(1/$3) lt 8 title "laenge=20",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaiv.txt" u (($2==36)?$1:1/0):(1/$3) lt 9 title "laenge=36",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaiv.txt" u (($2==50)?$1:1/0):(1/$3) lt 2 title "laenge=50",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaiv.txt" u (($2==100)?$1:1/0):(1/$3) lt 3 title "laenge=100",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaiv.txt" u (($2==200)?$1:1/0):(1/$3) lt 4 title "laenge=200",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaiv.txt" u (($2==350)?$1:1/0):(1/$3) lt 5 title "laenge=350",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaiv.txt" u (($2==500)?$1:1/0):(1/$3) lt 6 title "laenge=500"

set title "lcpunode02"
set yrange[0:4]
plot "Messungen/Zeiten/zeitenmehrerelaengenmittelnaivnode029cores.txt" u (($2==6)?$1:1/0):(1/$3):5 w yerrorbars lt 7 title "laenge=6",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaivnode029cores.txt" u (($2==10)?$1:1/0):(1/$3):5 w yerrorbars lt 1 title "laenge=10",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaivnode029cores.txt" u (($2==20)?$1:1/0):(1/$3):5 w yerrorbars lt 8 title "laenge=20",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaivnode029cores.txt" u (($2==36)?$1:1/0):(1/$3):5 w yerrorbars lt 9 title "laenge=36",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaivnode029cores.txt" u (($2==50)?$1:1/0):(1/$3):5 w yerrorbars lt 2 title "laenge=50",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaivnode029cores.txt" u (($2==100)?$1:1/0):(1/$3):5 w yerrorbars lt 3 title "laenge=100",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaivnode029cores.txt" u (($2==200)?$1:1/0):(1/$3):5 w yerrorbars lt 4 title "laenge=200",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaivnode029cores.txt" u (($2==350)?$1:1/0):(1/$3):5 w yerrorbars lt 5 title "laenge=350",\
	 "Messungen/Zeiten/zeitenmehrerelaengenmittelnaivnode029cores.txt" u (($2==500)?$1:1/0):(1/$3):5 w yerrorbars lt 6 title "laenge=500"
unset yrange
unset xrange

#set view map
#splot "Messungen/Zeiten/zeitenmessen-m001000-mehrerelaengenunddurchlaeufeqbig.txt" u 1:5:4 w image

set out "Messungen/speedup.pdf"

set title "Laenge=50"
plot "Messungen/Zeiten/zeitenmessen-laenge0050-m100000.txt" u 1:4 lt 7 title "" 
set title "Laenge=100"
plot "Messungen/Zeiten/zeitenmessen-laenge0100-m010000-mehrere.txt" u 1:4 lt 7 title "" 

set ylabel "Zeit(1 Core)/Zeit(n Cores)"

set title "Laenge=50"
plot "Messungen/Zeiten/zeitenmessen-laenge0050-m100000.txt" u 1:(1/$4) lt 7 title "" 
set title "Laenge=100"
plot "Messungen/Zeiten/zeitenmessen-laenge0100-m010000-mehrere.txt" u 1:(1/$4) lt 7 title "" 

set out "Messungen/speedup1000.pdf"
set title "Laenge=1000"
plot "Messungen/Zeiten/zeitenmessen-laenge1000-m001000-mehrere.txt" u 1:(1/$4) lt 7 title "" 

set ylabel "Zeit(n Cores)/Zeit(1 Core)"

set title "Laenge=1000"
plot "Messungen/Zeiten/zeitenmessen-laenge1000-m001000-mehrere.txt" u 1:4 lt 7 title "" 


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
