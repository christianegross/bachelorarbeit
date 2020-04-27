set ter pdfcairo size 4in,5.5in


#Magnetisierung
set out "Messungen/magnetisierung.pdf"
set xlabel "Temperatur"
set ylabel "Magnetisierung"

#set xrange [0:10.5]
#set yrange [-1.1:1.1]

#Ein Datenset
#plot "Messungen/Mittelwerte/messenmittelmagl0030.txt" u 1:2:3 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung"
#plot "Messungen/Mittelwerte/messenmittelmagl0300.txt" u 1:2:3 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung"

#Zwanzig Datensets
#set title "l=300 bis l=319"
#plot "Messungen/Mittelwerte/messenmittelmagl0301.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung"

#Zehn Datensets
#plot "Messungen/Mittelwerte/messenmittelmagl0200.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung"
plot "Messungen/Mittelwerte/messenmittel-l0050.txt" u 2:6:7 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung"

unset xrange
set xlabel "{/Symbol b}J"
plot "Messungen/Mittelwerte/messenmittel-l0050.txt" u 3:6:7 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung"
#plot "Messungen/Mittelwerte/messenmittelmagl0301.txt" u 3:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung"
#plot "Messungen/Mittelwerte/messenmittelmagl0200.txt" u 3:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung"

set xlabel "Temperatur"
set xrange [0:10.5]
#plot "Messungen/Mittelwerte/messenmittelmagl0100.txt" u 1:2:3 w yerrorbars ps 0.4 lt 100 title "laenge=100", "Messungen/Mittelwerte/messenmittelmagl0101.txt" u 1:2:3 w yerrorbars ps 0.4 lt 101 title "laenge=101", "Messungen/Mittelwerte/messenmittelmagl0102.txt" u 1:2:3 w yerrorbars ps 0.4 lt 102 title "laenge=102", "Messungen/Mittelwerte/messenmittelmagl0103.txt" u 1:2:3 w yerrorbars ps 0.4 lt 103 title "laenge=103", "Messungen/Mittelwerte/messenmittelmagl0104.txt" u 1:2:3 w yerrorbars ps 0.4 lt 104 title "laenge=104", "Messungen/Mittelwerte/messenmittelmagl0105.txt" u 1:2:3 w yerrorbars ps 0.4 lt 105 title "laenge=105", "Messungen/Mittelwerte/messenmittelmagl0106.txt" u 1:2:3 w yerrorbars ps 0.4 lt 106 title "laenge=106", "Messungen/Mittelwerte/messenmittelmagl0107.txt" u 1:2:3 w yerrorbars ps 0.4 lt 107 title "laenge=107", "Messungen/Mittelwerte/messenmittelmagl0108.txt" u 1:2:3 w yerrorbars ps 0.4 lt 108 title "laenge=108", "Messungen/Mittelwerte/messenmittelmagl0109.txt" u 1:2:3 w yerrorbars ps 0.4 lt 109 title "laenge=109"



#Akzeptanzrate
set out "Messungen/akzeptanzrate.pdf"
set xlabel "Temperatur"
set ylabel "Akzeptanzrate"
set yrange [-0.05:1.05]

#Ein Datenset
plot "Messungen/Mittelwerte/messenmittelakzl0300.txt" u 1:2:3 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate"
#Zwanzig Datensets
plot "Messungen/Mittelwerte/messenmittelakzl0301.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate"
plot "Messungen/Mittelwerte/messenmittel-l0050.txt" u 2:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Magnetisierung"
unset xrange
set xlabel "{/Symbol b}J"
#Zwanzig Datensets
plot "Messungen/Mittelwerte/messenmittelakzl0301.txt" u 3:4:5 w yerrorbars ps 0.4 lt 7 title "Messungen Akzeptanzrate"


set xlabel "Temperatur"
set xrange [0:10.5]

