set ter pdfcairo size 4in,5.5in


#Magnetisierung
set out "Messungen/magnetisierung.pdf"
set xlabel "Temperatur"
set ylabel "Magnetisierung"
set yrange [-0.05:1.05]

plot "Messungen/Mittelwerte/messenmittel-l0100.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=100"
plot "Messungen/Mittelwerte/messenmittel-l0010.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=10"
#plot "Messungen/Mittelwerte/messenmittel-l0050.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=50"
plot "Messungen/Mittelwerte/messenmittel-l0101.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=101"
plot "Messungen/Mittelwerte/messenmittel-l0102.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=102"
plot "Messungen/Mittelwerte/messenmittel-l0200.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=200"
#plot "Messungen/Mittelwerte/messenmittel-l0300.txt" u 2:6:7 w lines  lt 7 title "Messungen Magnetisierung, l=300"
plot "Messungen/Mittelwerte/messenmittel-l0101.txt" u 2:7 w lines  lt 7 title "Messfehler Magnetisierung, l=101"

set xlabel "{/Symbol b} J"

#plot "Messungen/Mittelwerte/messenmittel-l0100.txt" u 3:6:7 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung, l=100"
#plot "Messungen/Mittelwerte/messenmittel-l0101.txt" u 3:6:7 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung, l=101"
#plot "Messungen/Mittelwerte/messenmittel-l0200.txt" u 3:6:7 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung, l=200"
#plot "Messungen/Mittelwerte/messenmittel-l0300.txt" u 3:6:7 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung, l=300"




#Akzeptanzrate
set out "Messungen/akzeptanzrate.pdf"
set xlabel "Temperatur"
set ylabel "Akzeptanzrate"
set yrange [-0.05:1.05]

plot "Messungen/Mittelwerte/messenmittel-l0100.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=100"
plot "Messungen/Mittelwerte/messenmittel-l0100.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=10"
plot "Messungen/Mittelwerte/messenmittel-l0101.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=101"
plot "Messungen/Mittelwerte/messenmittel-l0102.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=102"
plot "Messungen/Mittelwerte/messenmittel-l0200.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=200"
plot "Messungen/Mittelwerte/messenmittel-l0300.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=300"


set xlabel "{/Symbol b} J"

plot "Messungen/Mittelwerte/messenmittel-l0100.txt" u 3:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=100"
plot "Messungen/Mittelwerte/messenmittel-l0101.txt" u 3:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=101"
plot "Messungen/Mittelwerte/messenmittel-l0200.txt" u 3:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=200"
plot "Messungen/Mittelwerte/messenmittel-l0300.txt" u 3:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate, l=300"

#Verschiedene arrays für vershiedene Längen, immer näher zusammen, um Verhalten vor allem bei niedrigen Temperaturen beobaachten zu können

