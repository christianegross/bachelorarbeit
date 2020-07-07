set ter epslatex size 6in, 4.5in color colortext

set style line 1 lc 7 lt 7 pt 7 #ps 0.2
set style line 2 lc 1 lt 7 pt 9
set style line 3 lc 2 lt 7 pt 11
set style line 4 lc 3 lt 7 pt 3
set style line 5 lc 4 lt 7 pt 5#4 gut

#Literaturergbnis Magnetisierung
litresult(x)=(1-(1/(sinh(2/x)*sinh(2/x)))**2)**(1.0/8)*(x<=2.269)+0*(x>2.269)

#Vergleich Hamiltonian
set out 'vergleichham.tex'
set xlabel 'Temperatur'
set ylabel '$H/\text{laenge}^2$'

set xrange [0:5]
set key top left
plot 'bootstrapalle-sweepaltohnepar-hamiltonian-l0120-m-010000-node00.txt' 			u (($2==128)?$6:1/0):($4/120**2):($5/120**2) w yerrorbars 	ls 4 title 'zeilenweise durchgehen',\
	 'bootstrapalle-sweepmehreregeneratoren-hamiltonian-l0120-m-010000-node00.txt' 	u (($2==128)?$6:1/0):($4/120**2):($5/120**2) w yerrorbars 	ls 3 title 'Schachbrettmuster parallel'
unset xrange
unset xlabel
unset ylabel

#Fehler beim bootstrapping, Temperatur=3
set out 'fehlerbootstrapping.tex'
set xlabel 'Blocklänge'
set ylabel 'Fehler'
set xrange[0.5:1540]
set logscale x

plot 'bootstrapfehler-magnetisierung-l0384-m-010000-node00.txt' u 2:5 ls 1 title '$\Delta M$',\
	 'bootstrapfehler-hamiltonian-l0384-m-010000-node00.txt' 	u 2:($5/384**2) ls 2 title '$\Delta H/\text{laenge^2}$',\
	 'bootstrapfehler-akzeptanz-l0384-m-010000-node00.txt' 		u 2:5 ls 3 title '$\Delta$ Akzeptanzrate',\

unset logscale x	 
unset xlabel
unset ylabel


set ter epslatex size 6in, 4.5in color colortext
#Vergleich Magnetisieurng verschiedene Längen
set xlabel 'Temperatur'
set ylabel '$M$'

set out 'magnetisierunglaenge.tex'
set key top right
set xrange [1:4]
set yrange[0:1]
plot 'bootstrapalle-magnetisierung-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 2 title '$\text{laenge}=120$',\
	 'bootstrapalle-magnetisierung-l0036-m-010000-node00.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 3 title '$\text{laenge}=36$',\
	 litresult(x) ls 1 title 'Theorie'
unset xrange
unset yrange

#ableitung
set out 'ableitung120128.tex'
set xlabel 'Temperatur'
set ylabel 'dM/dT'
set xrange[1:4]
plot 'ableitung-magnetisierung-laenge-0120-m-010000-node02.txt' u 1:2:3 w yerrorbars ls 1 title ''
unset xrange

#speedup
set out 'speeduplaenge.tex' 

einslinie(x)=1
set xlabel 'Verwendete Kerne'
set ylabel 'Speedup'
set key top left
set xrange[0:12.5]
plot 'zmittel-m001000-node02v5schedstatt05weninit-l10.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 1 title 'L=10',\
	 'zmittel-m001000-node02v5schedstatt05weninit-l40.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 2 title 'L=40',\
	 'zmittel-m001000-node02v5schedstatt05weninit-l70.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 3 title 'L=70',\
	 'zmittel-m001000-node02v5schedstatt05weninit-l500.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 4 title 'L=500',\
	 einslinie(x) lt 0 dt 5 title ''
unset yrange

set out 'speeduptemperatur.tex' 

einslinie(x)=1
set xlabel 'Verwendete Kerne'
set ylabel 'Speedup'
set key top left
set xrange[0:12.5]
plot 'zmittel-m001000-node02v5schedstatt05weninit-l500.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 4 title '$T=\num{0,5}$',\
	 'zmittel-m001000-node02v5schedstatt22weninit-l500.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 3 title '$T=\num{2,2}$',\
	 'zmittel-m001000-node02v5schedstatt23weninit-l500.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 2 title '$T=\num{2,3}$',\
	 'zmittel-m001000-node02v5schedstatt100weninit-l500.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 1 title '$T=\num{100}$',\
	 einslinie(x) lt 0 dt 5 title ''
unset yrange

set ter epslatex size 3.5in, 3.5in color colortext
#Schachbrett 

set out 'schachbrett.tex'
set view map
unset colorbox
unset xlabel
unset ylabel
set rmargin 0
set lmargin 0
set bmargin 0
set tmargin 0
set palette defined (0 'black' , 50 'gray', 100 'white')
set cbrange [0:1]
set xrange [-0.5:9.5] 
set yrange [9.5:-0.5] 
set title ''
splot 'schachbrett.txt' u 1:2:3 w image title ''

set output
