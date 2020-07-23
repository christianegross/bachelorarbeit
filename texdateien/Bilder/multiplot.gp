#Erstellung der Bilder fuer die bachelorarbeit, die als Multiplot ausgegeben werden, und deren nciht-Multiplot-Aequivalente
#Hier Multiplot mit 2x1Bildern von Magnetisierung, Akzeptanzrate und Hamiltonian, in grober und feiner Ansicht
#Zusaetzlich Plot von deren standardabweichungen

#Hamiltonian, Akzeptanzrate und Magnetisierung grob und fein
set ter epslatex size 3in, 3in color colortext

set macros
#In Anlehnung an http://www.gnuplotting.org/tag/multiplot/
LEFT='set lmargin at screen 0.15; set rmargin at screen 0.545; set xtics 0,1,4; set xrange [0:5];  set format y '
RIGHT='set lmargin at screen 0.555; set rmargin at screen 0.95; set xtics 0, 2000, 10000; set xrange [-10:10000]; unset ylabel; set format y '''

set style line 1 lc 7 lt 7 pt 7 #ps 0.2
set style line 2 lc 1 lt 7 pt 9
set style line 3 lc 2 lt 7 pt 11
set style line 4 lc 3 lt 7 pt 3
set style line 5 lc 4 lt 7 pt 5#4 gut

set xlabel 'Temperatur'

set out 'akzeptanzrategrob.tex'
set ylabel 'Akzeptanzrate'
set xtics 400
plot 'bootstrapalle-akzeptanz-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 1 title ''

set out 'magnetisierunggrob.tex'
set ylabel '$M$'
set xtics 400
plot 'bootstrapalle-magnetisierung-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 1 title ''

set out 'hamiltoniangrob.tex'
set ylabel '$H/\text{laenge}^2$'
set xtics 400
plot 'bootstrapalle-hamiltonian-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):($4/120**2):($5/120**2)  w yerrorbars ls 1 title ''
unset xtics
set xtics auto

set xrange[0:5]
set out 'akzeptanzratefein.tex'
set ylabel 'Akzeptanzrate'
plot 'bootstrapalle-akzeptanz-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 1 title ''

set out 'magnetisierungfein.tex'
set ylabel '$M$'
plot 'bootstrapalle-magnetisierung-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 1 title ''

set out 'hamiltonianfein.tex'
set ylabel '$H/\text{laenge}^2$'
plot 'bootstrapalle-hamiltonian-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):($4/120**2):($5/120**2)  w yerrorbars ls 1 title ''
unset xrange


set ter epslatex size 15 cm, 7.5cm color colortext
set out 'akzeptanzrate.tex'
set multiplot layout 1,2
set yrange [-0.01:1.01]
set ylabel 'Akzeptanzrate'

@LEFT
plot 'bootstrapalle-akzeptanz-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 1 title ''

@RIGHT
plot 'bootstrapalle-akzeptanz-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 1 title ''

unset multiplot

set out 'hamiltonian.tex'
set multiplot layout 1,2
set yrange [-2.01:0.01]
set ylabel '$H/\text{laenge}^2$'

@LEFT
plot 'bootstrapalle-hamiltonian-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):($4/120**2):($5/120**2)  w yerrorbars ls 1 title ''

@RIGHT
plot 'bootstrapalle-hamiltonian-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):($4/120**2):($5/120**2)  w yerrorbars ls 1 title ''

unset multiplot

set out 'magnetisierung.tex'
set multiplot layout 1,2
set yrange [-0.01:1.01]
set ylabel '$M$'

@LEFT
plot 'bootstrapalle-magnetisierung-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 1 title ''

@RIGHT
plot 'bootstrapalle-magnetisierung-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 1 title ''
unset multiplot

set output

set ter epslatex size 6in, 3in color colortext
set out 'akzeptanzratefehler.tex'
unset yrange
set multiplot layout 1,2
#set yrange [-0.01:1.01]
set ylabel 'Akzeptanzrate'

@LEFT
plot 'bootstrapalle-akzeptanz-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):5 ls 1 title ''

@RIGHT
plot 'bootstrapalle-akzeptanz-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):5  ls 1 title ''

unset multiplot

set out 'hamiltonianfehler.tex'
set multiplot layout 1,2
#set yrange [-2.01:0.01]
set ylabel '$H/\text{laenge}^2$'

@LEFT
plot 'bootstrapalle-hamiltonian-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):($5/120**2) ls 1 title ''

@RIGHT
plot 'bootstrapalle-hamiltonian-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):($5/120**2) ls 1 title ''

unset multiplot

set out 'hamiltonianfehler.tex'
set multiplot layout 1,2
#set yrange [-2.01:0.01]
set ylabel '$H/\text{laenge}^2$'

@LEFT
plot 'bootstrapalle-hamiltonian-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):($5/120**2) ls 1 title ''

@RIGHT
plot 'bootstrapalle-hamiltonian-l0120-m-010000-node02.txt' u (($2==128)?$6:1/0):($5/120**2) ls 1 title ''

unset multiplot

set out 'speeduptemperaturmultiplot.tex'
set multiplot layout 1,2
#set yrange [-0.01:1.01]
set ylabel 'Speedup'

@LEFT
einslinie(x)=1
ideal(x)=x
set xlabel 'OpenMP-Threads'
set ylabel 'Speedup'
set key top left
set xtics 0,2,10
set mxtics 2
set xrange[0:12.5]
set yrange[0:7]
plot 'zmittel-m001000-node02-l500.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 1 title '$T=\num{0,5}$',\
	 'zmittel-m001000-node02t22-l500.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 2 title '$T=\num{2,2}$',\
	 'zmittel-m001000-node02t23-l500.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 3 title '$T=\num{2,3}$',\
	 'zmittel-m001000-node02t100-l500.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 4 title '$T=\num{100}$',\
	 ideal(x) lt 0 dt 5 title ''#	 einslinie(x) lt 0 dt 5 title '',\
unset xrange
unset yrange

@RIGHT
set xlabel 'MPI-Prozesse'
set y2label 'Speedup'
set format y2 "$%h$"
set key top left
set xtics 0,5,20
set y2tics 0,5,20
set mxtics 5
set xrange[0:20.5]
set y2range[0:20.5]
set yrange[0:20.5]
plot 'skalqbigmpiverschtemp.txt' u (($4==0.5)?$3:1/0):10:11 w yerrorbars ls 1 axes x1y2 title '$T=\num{0,5}$',\
	 'skalqbigmpiverschtemp.txt' u (($4==2.2)?$3:1/0):10:11 w yerrorbars ls 2 axes x1y2 title '$T=\num{2,2}$',\
	 'skalqbigmpiverschtemp.txt' u (($4==2.3)?$3:1/0):10:11 w yerrorbars ls 3 axes x1y2 title '$T=\num{2,3}$',\
	 'skalqbigmpiverschtemp.txt' u (($4==100)?$3:1/0):10:11 w yerrorbars ls 4 axes x1y2 title '$T=\num{100}$',\
	 ideal(x) lt 0 dt 5 axes x1y2 title ''#einslinie(x) lt 0 dt 5 axes x1y2 title '',\
	 
unset xrange
unset y2range
unset multiplot

set output
