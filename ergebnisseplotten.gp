set ter pdfcairo size 4in,5.5in


#Magnetisierung
set out "Messungen/magnetisierung.pdf"
set xlabel "Temperatur"
set ylabel "Magnetisierung"
set yrange [-0.05:1.05]

#plot "Messungen/Mittelwerte/messenmittel-l0050-m-100000.txt" u 2:6:7 w yerrorbars ps 0.3  lt 7 title "l=50, 100000 Messungen"
plot "Messungen/Mittelwerte/messenmittel-l0050-m-001000.txt" u 2:6:7 w lines  lt 7 title "l=50, 1000 Messungen, neue Thermalisierung"
plot "Messungen/Mittelwerte/messenmittel-l0050-m-005000.txt" u 2:6:7 w lines  lt 7 title "l=50, 5000 Messungen, neue Thermalisierung"
set xrange [0.6:1.6]
plot "Messungen/Mittelwerte/messenmittel-l0080-m-001000.txt" u 2:6:7 w yerrorbars ps 0.3  lt 7 title "l=80, 1000 Messungen, neue Thermalisierung"
plot "Messungen/Mittelwerte/messenmittel-l0080-m-005000.txt" u 2:6:7 w yerrorbars ps 0.3  lt 7 title "l=80, 5000 Messungen, neue Thermalisierung"
unset xrange
#plot "Messungen/Mittelwerte/messenmittel-l0100-m-001000m.txt" u 2:6:7 w lines  lt 7 title "l=100, 1000 Messungen, mehrmalige Initialisierung"
#plot "Messungen/Mittelwerte/messenmittel-l0100-m-001000e.txt" u 2:6:7 w lines  lt 7 title "=100, 1000 Messungen, einmalige Initialisierung"
#plot "Messungen/Mittelwerte/messenmittel-l0100-m-010000m.txt" u 2:6:7 w lines  lt 7 title "l=100, 10000 Messungen, mehrmalige Initialisierung"
#plot "Messungen/Mittelwerte/messenmittel-l0010.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=10"
#plot "Messungen/Mittelwerte/messenmittel-l0050.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=50"
#plot "Messungen/Mittelwerte/messenmittel-l0051.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=51"
#plot "Messungen/Mittelwerte/messenmittel-l0005.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=5"
#plot "Messungen/Mittelwerte/messenmittel-l0101.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=101"
#plot "Messungen/Mittelwerte/messenmittel-l0102.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=102"
#plot "Messungen/Mittelwerte/messenmittel-l0200.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=200"
#plot "Messungen/Mittelwerte/messenmittel-l0300.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=300"
#plot "Messungen/Mittelwerte/messenmittel-l0101.txt" u 2:7 w lines  lt 7 title "Messfehler Magnetisierung, l=101"

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
set title "T=1,34, m=0,968885+/-0,134056"
#plot "Messungen/Messwerte/messung-l0050-t042.txt" u 1:3 ps 0.1 lt 1 title "Magnetisierung", "Messungen/Messwerte/messung-l0050-t042.txt" u 1:2 ps 0.1 lt 8 title "Akzeptanzrate",
set title "T=0,70, m=0.999978+/-0.000134"
#plot "Messungen/Messwerte/messung-l0050-t010.txt" u 1:3 ps 0.1 lt 1 title "Magnetisierung", "Messungen/Messwerte/messung-l0050-t042.txt" u 1:2 ps 0.1 lt 8 title "Akzeptanzrate",

#plot "Messungen/Messwerte/messung-l0080-t054.txt" u 1:3 ps 0.1 lt 1 title "Magnetisierung", "Messungen/Messwerte/messung-l0080-t054.txt" u 1:2 ps 0.1 lt 8 title "Akzeptanzrate",

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

