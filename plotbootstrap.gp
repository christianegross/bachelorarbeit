set ter pdfcairo size 4in,5.5in

set out "Messungen/bootstrap.pdf"

set logscale x
set xrange [1:1025]
set xlabel "l"

set title "T=0.015000"

set ylabel "Mittelwert"
f(x)=1.000000
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t000.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t000.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t000.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=0.165000"

set ylabel "Mittelwert"
f(x)=1.000000
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t010.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t010.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t010.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=0.315000"

set ylabel "Mittelwert"
f(x)=1.000000
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t020.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t020.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t020.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=0.465000"

set ylabel "Mittelwert"
f(x)=1.000000
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t030.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t030.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t030.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=0.615000"

set ylabel "Mittelwert"
f(x)=0.999995
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t040.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t040.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t040.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=0.765000"

set ylabel "Mittelwert"
f(x)=0.999940
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t050.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t050.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t050.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=0.915000"

set ylabel "Mittelwert"
f(x)=0.999665
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t060.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t060.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t060.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=1.065000"

set ylabel "Mittelwert"
f(x)=0.998789
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t070.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t070.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t070.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=1.215000"

set ylabel "Mittelwert"
f(x)=0.996741
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t080.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t080.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t080.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=1.365000"

set ylabel "Mittelwert"
f(x)=0.992721
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t090.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t090.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t090.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=1.515000"

set ylabel "Mittelwert"
f(x)=0.985571
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t100.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t100.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t100.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=1.665000"

set ylabel "Mittelwert"
f(x)=0.973743
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t110.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t110.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t110.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=1.815000"

set ylabel "Mittelwert"
f(x)=0.954374
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t120.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t120.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t120.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=1.965000"

set ylabel "Mittelwert"
f(x)=0.921954
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t130.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t130.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t130.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=2.115000"

set ylabel "Mittelwert"
f(x)=0.859346
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t140.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t140.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t140.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=2.265000"

set ylabel "Mittelwert"
f(x)=0.645634
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t150.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t150.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t150.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=2.415000"

set ylabel "Mittelwert"
f(x)=0.191903
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t160.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t160.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t160.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=2.565000"

set ylabel "Mittelwert"
f(x)=0.101350
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t170.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t170.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t170.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=2.715000"

set ylabel "Mittelwert"
f(x)=0.077442
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t180.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t180.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t180.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=2.865000"

set ylabel "Mittelwert"
f(x)=0.062966
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t190.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t190.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t190.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=3.015000"

set ylabel "Mittelwert"
f(x)=0.054531
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t200.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t200.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t200.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=3.165000"

set ylabel "Mittelwert"
f(x)=0.047957
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t210.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t210.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t210.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=3.315000"

set ylabel "Mittelwert"
f(x)=0.042772
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t220.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t220.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t220.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=3.465000"

set ylabel "Mittelwert"
f(x)=0.040025
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t230.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t230.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t230.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=3.615000"

set ylabel "Mittelwert"
f(x)=0.037902
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t240.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t240.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t240.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=3.765000"

set ylabel "Mittelwert"
f(x)=0.035642
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t250.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t250.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t250.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=3.915000"

set ylabel "Mittelwert"
f(x)=0.033973
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t260.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t260.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t260.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=4.065000"

set ylabel "Mittelwert"
f(x)=0.032619
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t270.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t270.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t270.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=4.215000"

set ylabel "Mittelwert"
f(x)=0.031119
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t280.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t280.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t280.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"

set title "T=4.365000"

set ylabel "Mittelwert"
f(x)=0.030459
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t290.txt" u 1:3:4 w yerrorbars lt 7 ps 0.3 title "Bootstrap-Mittelwerte", f(x) lt 6 title "naiver Mittelwert"

set ylabel "Varianz"
plot "Messungen/Bootstrapwerte/bootstrap-laenge0050-t290.txt" u 1:4 lt 7 ps 0.4 title "Bootstrap-Varianz", "Messungen/Bootstrapwerte/bootstrap-laenge0050-t290.txt" u 1:4 w lines lt 7 title "Bootstrap-Varianz"
