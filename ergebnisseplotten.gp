set ter pdfcairo size 4in,5.5in


#Magnetisierung
set out "Messungen/magnetisierung.pdf"
set xlabel "Temperatur"
set ylabel "Magnetisierung"
set yrange [-0.05:1.05]

#plot "Messungen/Mittelwerte/messenmittel-l0050-m-100000.txt" u 2:6:7 w yerrorbars ps 0.3  lt 7 title "l=50, 100000 Messungen"
#plot "Messungen/Mittelwerte/messenmittel-l0050-m-001000.txt" u 2:6:7 w lines  lt 7 title "l=50, 1000 Messungen, neue Thermalisierung"


plot "Messungen/Mittelwerte/messenmittel-l0050-m-005000.txt" u 2:6:7 w lines  lt 7 title "l=50, 5000 Messungen, neue Thermalisierung"
plot "Messungen/Mittelwerte/messenmittel-l0050-m-005000.txt" u 2:6:7 w yerrorbars ps 0.3  lt 7 title "l=50, 5000 Messungen, neue Thermalisierung"
plot "Messungen/Mittelwerte/messenmittel-l0080-m-005000.txt" u 2:6:7 w yerrorbars ps 0.3  lt 7 title "l=80, 5000 Messungen, neue Thermalisierung"
plot "Messungen/Mittelwerte/messenmittel-l0080-m-005000.txt" u 2:6:7 w lines  lt 7 title "l=80, 5000 Messungen, neue Thermalisierung"

plot "Messungen/Mittelwerte/messenmittel-l0010-m-005000.txt" u 2:6:7 w yerrorbars ps 0.3  lt 7 title "l=10, 5000 Messungen, neue Thermalisierung"
plot "Messungen/Mittelwerte/messenmittel-l0010-m-005000.txt" u 2:6:7 w lines  lt 7 title "l=10, 5000 Messungen, neue Thermalisierung"
unset xrange

#set xrange [0.6:1.6]
#plot "Messungen/Mittelwerte/messenmittel-l0080-m-001000.txt" u 2:6:7 w yerrorbars ps 0.3  lt 7 title "l=80, 1000 Messungen, neue Thermalisierung"
#unset xrange

set xlabel "{/Symbol b} J"

#plot "Messungen/Mittelwerte/messenmittel-l0100.txt" u 3:6:7 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung, l=100"
#plot "Messungen/Mittelwerte/messenmittel-l0101.txt" u 3:6:7 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung, l=101"
#plot "Messungen/Mittelwerte/messenmittel-l0200.txt" u 3:6:7 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung, l=200"
#plot "Messungen/Mittelwerte/messenmittel-l0300.txt" u 3:6:7 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung, l=300"


set out "thermalisierung.pdf"
#set xrange [0:10000]
set xlabel "Messung"
set ylabel "Messwert"

set key center right
#set title "T=1,34, m=0,968885+/-0,134056"
plot "Messungen/Messwerte/messung-l0050-t150.txt" u 1:3 ps 0.1 lt 1 title "Magnetisierung", "Messungen/Messwerte/messung-l0050-t150.txt" u 1:2 ps 0.1 lt 8 title "Akzeptanzrate",
#set title "T=0,70, m=0.999978+/-0.000134"
#plot "Messungen/Messwerte/messung-l0050-t010.txt" u 1:3 ps 0.1 lt 1 title "Magnetisierung", "Messungen/Messwerte/messung-l0050-t042.txt" u 1:2 ps 0.1 lt 8 title "Akzeptanzrate",

plot "Messungen/Messwerte/messung-l0080-t152.txt" u 1:3 ps 0.1 lt 1 title "Magnetisierung", "Messungen/Messwerte/messung-l0080-t152.txt" u 1:2 ps 0.1 lt 8 title "Akzeptanzrate",
plot "Messungen/Messwerte/messung-l0080-t145.txt" u 1:3 ps 0.1 lt 1 title "Magnetisierung", "Messungen/Messwerte/messung-l0080-t145.txt" u 1:2 ps 0.1 lt 8 title "Akzeptanzrate",

unset xrange

#Akzeptanzrate
set out "Messungen/akzeptanzrate.pdf"
set xlabel "Temperatur"
set ylabel "Akzeptanzrate"
set yrange [-0.05:1.05]

plot "Messungen/Mittelwerte/messenmittel-l0050-m-005000.txt" u 2:4:5 w yerrorbars ps 0.3  lt 7 title "l=50, 5000 Messungen"
#plot "Messungen/Mittelwerte/messenmittel-l0100.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=100"
#plot "Messungen/Mittelwerte/messenmittel-l0100.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=10"
#plot "Messungen/Mittelwerte/messenmittel-l0101.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=101"
#plot "Messungen/Mittelwerte/messenmittel-l0102.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=102"
#plot "Messungen/Mittelwerte/messenmittel-l0200.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=200"
#plot "Messungen/Mittelwerte/messenmittel-l0300.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=300"


set xlabel "{/Symbol b} J"

#plot "Messungen/Mittelwerte/messenmittel-l0100.txt" u 3:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=100"
#plot "Messungen/Mittelwerte/messenmittel-l0101.txt" u 3:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=101"
#plot "Messungen/Mittelwerte/messenmittel-l0200.txt" u 3:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=200"
#plot "Messungen/Mittelwerte/messenmittel-l0300.txt" u 3:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=300"

#Verschiedene arrays für vershiedene Längen, immer näher zusammen, um Verhalten vor allem bei niedrigen Temperaturen beobaachten zu können

