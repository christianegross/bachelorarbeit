#plot der Bilder, die fuer Bachelorarbeit benoetigt werden und nicht multiplot sind, z.B. Magnetisierung bei verschiedenen Laengen, Hamiltonian von verschiedenen sweep-Funktionen, speedup von cores, tkrit(l), schachbrettmuster
set ter epslatex size 15 cm, 10.6cm color colortext

set style line 1 lc 7 lt 7 pt 7 #ps 0.2
set style line 2 lc 1 lt 7 pt 8
set style line 3 lc 2 lt 7 pt 11
set style line 4 lc 3 lt 7 pt 3
set style line 5 lc 4 lt 7 pt 5#4 gut
set style line 6 lc 5 lt 7 pt 4#4 gut


#Literaturergbnis Magnetisierung
litresult(x)=(1-(1/(sinh(2/x)*sinh(2/x)))**2)**(1.0/8)*(x<=2.269)+0*(x>2.269)

#Vergleich Observablen
set out 'vergleichham.tex'
set xlabel 'Temperatur'
set ylabel '$H/\text{laenge}^2$'

set xrange [0:5]
set key top left
plot 'keinpar-bootstrapalle-hamiltonian-l0100-m-010240-node00-sch-10.txt' every 10 			u (($2==128)?$6:1/0):($4/100**2):($5/100**2) w yerrorbars 	ls 4 title 'seriell',\
	 'bootstrapalle-hamiltonian-l0100-m-010240-node02-sch-01.txt' 			every 100		u (($2==128)?$6:1/0):($4/100**2):($5/100**2) w yerrorbars 	ls 3 title 'OpenMP',\
	 'bootstrapalle-hamiltonian-l0100-m-010000-proz10-sch-01.txt' 			every 100		u (($2==128)?$6:1/0):($4/100**2):($5/100**2) w yerrorbars 	ls 2 title 'MPI'
unset xrange
unset xlabel
unset ylabel

set out 'vergleichmag.tex'
set xlabel 'Temperatur'
set ylabel '$M$'

set xrange [0:5]
set key top left
plot 'keinpar-bootstrapalle-magnetisierung-l0100-m-010240-node00-sch-10.txt' every 10 			u (($2==128)?$6:1/0):($4/100**2):($5/100**2) w yerrorbars 	ls 4 title 'seriell',\
	 'bootstrapalle-magnetisierung-l0100-m-010240-node02-sch-01.txt' 			every 100		u (($2==128)?$6:1/0):($4/100**2):($5/100**2) w yerrorbars 	ls 3 title 'OpenMP',\
	 'bootstrapalle-magnetisierung-l0100-m-010000-proz10-sch-01.txt' 			every 100		u (($2==128)?$6:1/0):($4/100**2):($5/100**2) w yerrorbars 	ls 2 title 'MPI'
unset xrange
unset xlabel
unset ylabel

set out 'vergleichakz.tex'
set xlabel 'Temperatur'
set ylabel 'Akzeptanzrate'

set xrange [0:5]
set key top left
plot 'keinpar-bootstrapalle-akzeptanz-l0100-m-010240-node00-sch-10.txt' every 10 			u (($2==128)?$6:1/0):($4/100**2):($5/100**2) w yerrorbars 	ls 4 title 'seriell',\
	 'bootstrapalle-akzeptanz-l0100-m-010240-node02-sch-01.txt' 			every 100		u (($2==128)?$6:1/0):($4/100**2):($5/100**2) w yerrorbars 	ls 3 title 'OpenMP',\
	 'bootstrapalle-akzeptanz-l0100-m-010000-proz10-sch-01.txt' 			every 100		u (($2==128)?$6:1/0):($4/100**2):($5/100**2) w yerrorbars 	ls 2 title 'MPI'
unset xrange
unset xlabel
unset ylabel


set ter epslatex size 15 cm, 10.6cm color colortext
#Vergleich Magnetisieurng verschiedene Längen
set xlabel 'Temperatur'
set ylabel '$M$'

set out 'magnetisierunglaenge.tex'
set key top right
set xrange [1:3.6]
set yrange[0:1]
plot litresult(x) ls 1 lw 2 title 'Theorie',\
	 'bootstrapalle-magnetisierung-l0012-m-010240-node00-sch-01.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 2 title '$L=12$',\
	 'bootstrapalle-magnetisierung-l0024-m-010240-node00-sch-01.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 3 title '$L=24$',\
	 'bootstrapalle-magnetisierung-l0048-m-010240-node00-sch-01.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 4 title '$L=48$',\
	 'bootstrapalle-magnetisierung-l0768-m-010240-node00-sch-01.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 5 title '$L=768$'
unset xrange
unset yrange

set out 'magnetisierunglaengempi.tex'
set key top right
set xrange [1:3.6]
set yrange[0:1]
plot litresult(x) ls 1 lw 2 title 'Theorie',\
	 'bootstrapalle-magnetisierung-l0012-m-010000-proz02-sch-01.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 2 title '$L=12$',\
	 'bootstrapalle-magnetisierung-l0024-m-010000-proz02-sch-01.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 3 title '$L=24$',\
	 'bootstrapalle-magnetisierung-l0048-m-010000-proz06-sch-01.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 4 title '$L=48$',\
	 'bootstrapalle-magnetisierung-l0480-m-010000-proz20-sch-01.txt' u (($2==128)?$6:1/0):4:5  w yerrorbars ls 5 title '$L=480$'
unset xrange
unset yrange


#speedup
set out 'speeduplaenge.tex' 

einslinie(x)=1
set xlabel 'Verwendete Threads'
set ylabel 'Speedup'
set key top left
set mxtics 2
set xrange[0:12.5]
plot 'zmittel-m001000-node02-l10.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 1 title 'L=10',\
	 'zmittel-m001000-node02-l40.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 2 title 'L=40',\
	 'zmittel-m001000-node02-l110.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 3 title 'L=110',\
	 'zmittel-m001000-node02-l500.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 4 title 'L=500',\
	 'zmittel-m001000-node00skalopenmp-l1260.txt' u (($1!=0)?$1:1/0):(1/$5):($6/$5/$5) w yerrorbars ls 5 title 'L=1260',\
	 einslinie(x) lt 0 dt 1 title ''
unset yrange
unset mxtics


set out 'speeduplaengempi.tex' 

einslinie(x)=1
set xlabel 'Anzahl MPI-Prozesse'
set ylabel 'Speedup'
set key top left
set xrange[0:20.5]
set mxtics 5
plot 'skalqbigmpiverschlaengen.txt' u (($5==12)?$3:1/0):10:11 w yerrorbars  ls 1 title 'L=12',\
	 'skalqbigmpiverschlaengen.txt' u (($5==60)?$3:1/0):10:11 w yerrorbars  ls 2 title 'L=60',\
	 'skalqbigmpiverschlaengen.txt' u (($5==240)?$3:1/0):10:11 w yerrorbars  ls 3 title 'L=240',\
	 'skalqbigmpiverschlaengen.txt' u (($5==840)?$3:1/0):10:11 w yerrorbars  ls 4 title 'L=840',\
	 'skalqbigmpiverschlaengen.txt' u (($5==1260)?$3:1/0):10:11 w yerrorbars  ls 5 title 'L=1260',\
	 einslinie(x) lt 0 dt 5 title ''
unset xrange
unset mxtics


set out 'tkritvonl.tex'
tkrit(x)=2.269
yfehler=0.002
set logscale x
set key top right
set xrange [10:1000]
set xtics 12, 2, 1000
set xlabel 'Gitterlänge'
set ylabel '$T_c$'
plot tkrit(x) lt 0 dt 5 title 'theoretischer Wert',\
	 'tkritopenmp25l128.txt' u ((($1==1)&&($3==2.0)&&($4==0))?$2:1/0):6:(0.002) w yerrorbars ls 1 title 'bestimmte Werte OpenMP',\
	 'tkritmpi15l128.txt' u ((($1==1)&&($3==2.0)&&($4==0))?$2:1/0):6:(0.002) w yerrorbars ls 2 title 'bestimmte Werte MPI'
unset logscale x
unset xrange
unset mxtics

set ter epslatex size 3in, 3in color colortext
#Schachbrett 

set out 'schachbrett.tex'
set view map
unset colorbox
unset xlabel
unset ylabel
set ytics out nomirror
set xtics out nomirror 0,1,8
set rmargin 0
set lmargin 0
set bmargin 0
set tmargin 0
set palette defined (0 'black' , 50 'gray', 100 'white')
set cbrange [0:1]
set xrange [-0.5:7.5] 
set yrange [7.5:-0.5] 
set title ''
splot 'schachbrett.txt' u 1:2:3 w image title ''

set output
