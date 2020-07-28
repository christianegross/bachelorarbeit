set ter pdfcairo size 4in,5.5in

set out "Messungen/magmin0407.pdf"

set xlabel "L"
set ylabel "Tkrit"
tkrit(x)=2.269
einhuellende(x)=2.269+a0+a1*exp(-1*a2*x)
a0=0.01
a1=0.8
a2=0.25
fit einhuellende(x) "Messungen/Ableitungen/mins/25laengenl128.txt" u ((($1==1)&&($3==2.0)&&($4==0))?$2:1/0):6 via a0, a1, a2

unset xrange
set logscale x
plot tkrit(x) dt 2 lc 0 title "theorie",\
	 "Messungen/Ableitungen/mins/25laengenl128.txt" u ((($1==1)&&($3==2.0)&&($4==0))?$2:1/0):6 lt 7 title "2, norm, 128",\
	 "Messungen/Ableitungen/mins/25laengenl256.txt" u ((($1==1)&&($3==2.0)&&($4==0))?$2:1/0):6 lt 8 title "2, norm, 256",\
	 einhuellende(x) title "einhuellende"
plot tkrit(x) dt 2 lc 0 title "theorie",\
	 "Messungen/Ableitungen/mins/25laengenl128.txt" u ((($1==1)&&($3==2.0)&&($4==0))?$2:1/0):6 lt 7 title "2, norm",\
	 "Messungen/Ableitungen/mins/25laengenl128.txt" u ((($1==1)&&($3==3.0)&&($4==0))?$2:1/0):6 lt 8 title "3, norm",\
	 "Messungen/Ableitungen/mins/25laengenl128.txt" u ((($1==1)&&($3==2.0)&&($4==1))?$2:1/0):6 lt 9 title "2, sqrt",\
	 "Messungen/Ableitungen/mins/25laengenl128.txt" u ((($1==1)&&($3==3.0)&&($4==1))?$2:1/0):6 lt 10 title "3, sqrt"
plot tkrit(x) dt 2 lc 0 title "theorie",\
	 "Messungen/Ableitungen/mins/allesm10368mit768mitnormundsqrt.txt" u ((($1==1)&&($3==2.0)&&($4==0))?$2:1/0):6 lt 7 title "2, norm",\
	 "Messungen/Ableitungen/mins/allesm10368mit768mitnormundsqrt.txt" u ((($1==1)&&($3==3.0)&&($4==0))?$2:1/0):6 lt 8 title "3, norm",\
	 "Messungen/Ableitungen/mins/allesm10368mit768mitnormundsqrt.txt" u ((($1==1)&&($3==2.0)&&($4==1))?$2:1/0):6 lt 9 title "2, sqrt",\
	 "Messungen/Ableitungen/mins/allesm10368mit768mitnormundsqrt.txt" u ((($1==1)&&($3==3.0)&&($4==1))?$2:1/0):6 lt 10 title "3, sqrt"
unset logscale x
#plot tkrit(x) dt 2 lc 0 title "theorie",\
#	 "Messungen/Ableitungen/mins/gesm10368norm.txt" u ((($1==1)&&($3==2.0))?$2:1/0):5 lt 7 title "sch 1, 2",\
#	 "Messungen/Ableitungen/mins/gesm10368norm.txt" u ((($1==1)&&($3==3.0))?$2:1/0):5 lt 8 title "sch 1, 3",\
#	 "Messungen/Ableitungen/mins/gesm10368norm.txt" u ((($1==2)&&($3==2))?$2:1/0):5 lt 9 title "sch 2, 2",\
#	 "Messungen/Ableitungen/mins/gesm10368norm.txt" u ((($1==2)&&($3==3))?$2:1/0):5 lt 10 title "sch 2, 3",\
#	 "Messungen/Ableitungen/mins/gesm10368norm.txt" u ((($1==5)&&($3==2))?$2:1/0):5 lt 11 title "sch 5, 2",\
#	 "Messungen/Ableitungen/mins/gesm10368norm.txt" u ((($1==5)&&($3==3))?$2:1/0):5 lt 12 title "sch 5, 3",\
#	 "Messungen/Ableitungen/mins/gesm10368norm.txt" u ((($1==10)&&($3==2))?$2:1/0):5 lt 13 title "sch 10, 2",\
#	 "Messungen/Ableitungen/mins/gesm10368norm.txt" u ((($1==10)&&($3==3))?$2:1/0):5 lt 14 title "sch 10, 3"
	 
#plot tkrit(x) dt 2 lc 0 title "theorie",\
#	 "Messungen/Ableitungen/mins/gesm10368sqrt.txt" u ((($1==1)&&($3==2))?$2:1/0):5 lt 7 title "sch 1, 2",\
#	 "Messungen/Ableitungen/mins/gesm10368sqrt.txt" u ((($1==1)&&($3==3))?$2:1/0):5 lt 8 title "sch 1, 3",\
#	 "Messungen/Ableitungen/mins/gesm10368sqrt.txt" u ((($1==2)&&($3==2))?$2:1/0):5 lt 9 title "sch 2, 2",\
#	 "Messungen/Ableitungen/mins/gesm10368sqrt.txt" u ((($1==2)&&($3==3))?$2:1/0):5 lt 10 title "sch 2, 3",\
#	 "Messungen/Ableitungen/mins/gesm10368sqrt.txt" u ((($1==5)&&($3==2))?$2:1/0):5 lt 11 title "sch 5, 2",\
#	 "Messungen/Ableitungen/mins/gesm10368sqrt.txt" u ((($1==5)&&($3==3))?$2:1/0):5 lt 12 title "sch 5, 3",\
#	 "Messungen/Ableitungen/mins/gesm10368sqrt.txt" u ((($1==10)&&($3==2))?$2:1/0):5 lt 13 title "sch 10, 2",\
#	 "Messungen/Ableitungen/mins/gesm10368sqrt.txt" u ((($1==10)&&($3==3))?$2:1/0):5 lt 14 title "sch 10, 3"
unset xlabel
unset ylabel

set out "Messungen/mag0107.pdf"

#set xrange [2.2:3]
plot "testtemp.txt" u 0:1 w lp lt 7 ps 0.3
unset xrange
plot "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0012-m-010240-node02.txt" using (($2==128)?$6:1/0):4:5 w yerrorbars lt 7 ps 0.3 title "12",\
  "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0024-m-010240-node02.txt" using (($2==128)?$6:1/0):4:5 w yerrorbars lt 8 ps 0.3 title "24",\
  "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0036-m-010240-node02.txt" using (($2==128)?$6:1/0):4:5 w yerrorbars lt 11 ps 0.3 title "36",\
  "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0048-m-010240-node02.txt" using (($2==128)?$6:1/0):4:5 w yerrorbars lt 12 ps 0.3 title "48",\
  "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0096-m-010240-node02.txt" using (($2==128)?$6:1/0):4:5 w yerrorbars lt 9 ps 0.3 title "96",\
  "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0192-m-010240-node02.txt" using (($2==128)?$6:1/0):4:5 w yerrorbars lt 10 ps 0.3 title "192"
  
plot "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0192-m-010240-node02.txt" using (($2==128)?$6:1/0):4 lt 10 ps 0.3 title "192"
  
plot "Messungen/ableitungdreipunkt-magnetisierung-laenge-0012-m-010240-node02.txt" using 1:2 lt 7 ps 0.3 title "l=12"
plot "Messungen/ableitungdreipunkt-magnetisierung-laenge-0192-m-010240-node02.txt" using 1:2 lt 10 ps 0.3 title "l=192"

plot "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0040-m-010240-node00.txt" using (($2==128)?$6:1/0):4:5 w yerrorbars lt 7 ps 0.3 title ""#, litresult(x) lt 6 title ""
plot "Messungen/Bootstrapges/bootstrapalle-magquad-l0040-m-010240-node00.txt" using (($2==128)?$6:1/0):(sqrt($4)):($5/sqrt($4)) w yerrorbars lt 7 ps 0.3 title ""#, litresult(x) lt 6 title ""
plot "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0070-m-020480-node02.txt" using (($2==128)?$6:1/0):4:5 w yerrorbars lt 7 ps 0.3 title ""#, litresult(x) lt 6 title ""

litresult(x)=(1-(1/(sinh(2/x)*sinh(2/x)))**2)**(1.0/8)
set out "Messungen/wurzelmagquad.pdf"
plot "Messungen/Bootstrapges/bootstrapalle-magquad-l0048-m-010000-node00.txt" u (($2==128)?$6:1/0):($4)**(0.5):($5)/(($4)**(0.5)) lt 7 ps 0.4 title "1", "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0048-m-010000-node00.txt" u (($2==128)?$6:1/0):4:5 lt 1 ps 0.2 title "3" , litresult(x) lt 6 title "2"
#sqrt($4):($5/sqrt($4))

set out "Messungen/ergebnisseparallel.pdf"
set xlabel "Temperatur"
set ylabel "Magnetisierung"
set yrange [-0.05:1.05]
set xrange [0:7.5]
set title "laenge=120"
plot "Messungen/Bootstrapges/bootstrapalle-l0120-m-010000-node02.txt" using (($2==128)?$6:1/0):4:5 w yerrorbars lt 7 ps 0.3 title "", litresult(x) lt 6 title ""
unset yrange
#(($1==128)?$5:1/0)
set ylabel "Ableitung Magnetisierung"
plot "Messungen/Ableitungen/ableitung-laenge-0120-m-010000-node02.txt" u 1:2:3 w yerrorbars lt 7 ps 0.3 title ""

set yrange [-0.05:1.05]
set xrange [1.5:3]
set title "laenge=120"
plot "Messungen/Bootstrapges/bootstrapalle-l0120-m-010000-node02.txt" using (($2==128)?$6:1/0):4:5 w yerrorbars lt 7 ps 0.3 title "", litresult(x) lt 6 title ""
unset yrange
set ylabel "Ableitung Magnetisierung"
plot "Messungen/Ableitungen/ableitung-laenge-0120-m-010000-node02.txt" u 1:2:3 w yerrorbars lt 7 ps 0.3 title ""

set title "Messungen QBiG"
unset xrange
set ylabel "Akzeptanzrate"
plot "Messungen/Bootstrapges/bootstrapalle-akzeptanz-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5  w yerrorbars ps 0.3 lt 7 title ""
set xrange [0:5]
set ylabel "Akzeptanzrate"
plot "Messungen/Bootstrapges/bootstrapalle-akzeptanz-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5  w yerrorbars ps 0.3 lt 7 title ""
set xrange [05:20]
set ylabel "Akzeptanzrate"
plot "Messungen/Bootstrapges/bootstrapalle-akzeptanz-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5  w yerrorbars ps 0.3 lt 7 title ""

unset xrange
set ylabel "hamiltonian*laenge^2"
plot "Messungen/Bootstrapges/bootstrapalle-hamiltonian-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):($4/120**2):($5/120**2) w yerrorbars ps 0.3 lt 7 title ""
set xrange [0:5]
set ylabel "hamiltonian*laenge^2"
plot "Messungen/Bootstrapges/bootstrapalle-hamiltonian-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):($4/120**2):($5/120**2) w yerrorbars ps 0.3 lt 7 title ""
set xrange [05:20]
set ylabel "hamiltonian*laenge^2"
plot "Messungen/Bootstrapges/bootstrapalle-hamiltonian-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):($4/120**2):($5/120**2) w yerrorbars ps 0.3 lt 7 title ""

unset xrange
set ylabel "magnetisierung"
plot "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5 w yerrorbars ps 0.3 lt 7 title ""
set xrange [0:5]
set ylabel "magnetisierung"
plot "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5 w yerrorbars ps 0.3 lt 7 title ""
set xrange [05:20]
set ylabel "magnetisierung"
plot "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5 w yerrorbars ps 0.3 lt 7 title ""

unset xrange
set ylabel "rel. Fehler Akzeptanzrate"
plot "Messungen/Bootstrapges/bootstrapalle-akzeptanz-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):($5/$4) ps 0.3 lt 7 title ""

set ylabel "rel. Fehler Magnetisierung"
plot "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):($5/$4) ps 0.3 lt 7 title ""

set ylabel "rel. Fehler Hamiltonian"
plot "Messungen/Bootstrapges/bootstrapalle-hamiltonian-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):($5/$4) ps 0.3 lt 7 title ""


set xrange [0:5]
set ylabel "magnetisierung"
plot "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0120-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5 w yerrorbars ps 0.3 lt 7 title "120",\
 "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0048-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5 w yerrorbars ps 0.3 lt 6 title "48",\
 "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0480-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5 w yerrorbars ps 0.3 lt 6 title "480",\
 "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0002-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5 w yerrorbars ps 0.3 lt 4 title "2",\
 "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0012-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5 w yerrorbars ps 0.3 lt 5 title "12"
 
 
plot "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0480-m-010000-node02.txt" u 6:4 w lp lt 6 title "480", \
  "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0480-m-010000-node02.txt" u (($2==128)?$6:1/0):4:5 w yerrorbars ps 0.3 lt 6 title "480"
