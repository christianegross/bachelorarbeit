set terminal epslatex size 6in, 4.5in color colortext

set style line 1 lc 1 pt 4
set style line 2 lc 2 pt 8 
set style line 3 lc 3 pt 10

#Vergleichsweepparallelisierung
set output "vergleichham.tex"
set title "Hamiltonian"
set xlabel "Temperatur"
set ylabel "H"
set key bottom right

dateimittel="vergleichmittel-l0050-m-010000.txt"

plot dateimittel u 1:2:3 w yerrorbars ls 1   title "Zeilenweise",\
dateimittel u 1:4:5 w yerrorbars ls 2   title "Schachbrett nicht parallel",\
dateimittel u 1:6:7 w yerrorbars ls 3   title "Schachbrett parallel",\


set output "vergleichakz.tex"
set title "Akzeptanzrate"
set ylabel "Akzeptanzrate"

plot dateimittel u 1:8:9 w yerrorbars ls 1   title "Zeilenweise",\
dateimittel u 1:10:11 w yerrorbars ls 2   title "Schachbrett nicht parallel",\
dateimittel u 1:12:13 w yerrorbars ls 3   title "Schachbrett parallel",\

#naive Ergebnisse
set xlabel "Temperatur"
set ylabel "Magnetisierung"
set yrange [-0.05:1.05]
set key top right


set out "magnetisierungnaiv.tex"
set title "l=50"
plot "messenmittel-l0050-m-010000.txt" u 2:6:7 w yerrorbars lt 7 title ""

set out "akzeptanzratenaiv.tex"
plot "messenmittel-l0050-m-010000.txt" u 2:4:5 w yerrorbars  lt 7 title ""
unset title

#MagnetisierungbeiBootstrap, Bootstrapmitundohnpar

set out "magnetisierungbootstrap.tex"
do for [n in "32 64 128 256 384 512 640 758 876 1024 1280 1536"]{
ausgabe=sprintf("magnetisierungbootstrap-l-%.4s.tex", n)
set out ausgabe
titel ="l=".n
plot "bootstrapalle-l0050-m-010000-nurpar.txt" using ((($2==n)&&($1==1))?$6:1/0):4:5 w yerrorbars lt 7 title titel
}


do for [n in "32 64 128 256 384 512 640 758 876 1024 1280 1536"]{
ausgabe=sprintf("vergleichbootstrap-l-%.4s.tex", n)
set out ausgabe
set title "l=".n
plot "bootstrapalle-l0050-m-010000.txt" using ((($2==n)&&($1==1))?$6:1/0):4:5 w yerrorbars lt 7 title "mit par",\
	"bootstrapalle-l0050-m-010000.txt" using ((($2==n)&&($1==2))?$6:1/0):4:5 w yerrorbars pt 10 lc 6 title "ohne par"
}

#Gitter
set view map
unset xlabel
unset ylabel
set ter epslatex size 4in, 4in color colortext

set palette defined (0 "white" , 50 "gray", 100 "black")
set title "Thermalisiert" font "arial,40"
set cbrange[-1:1]
unset colorbox

set xrange [-0.5:49.5]
set yrange [-0.5:49.5]
do for [t=0:295:5]{
ausgabe=sprintf("Gitter/gitter-laenge0050-t%.3d.tex", t)
datei=sprintf("Gitter/thermalisierung-laenge0050-t%.3d.txt", t)
set out ausgabe
titleplot=sprintf("T=%f", 0.015+t*0.015)
set title titleplot
splot datei u 1:2:3 w image title "" 
}
