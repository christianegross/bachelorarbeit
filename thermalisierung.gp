set ter pdfcairo size 8in, 11in


set out "thermalisierungplot.pdf"
plot "thermalisierung.txt" u 1:2 lt 7 w lines


set out "messungakzeptanzrate.pdf"
set yrange [0:1]
set xlabel "Messung"
set ylabel "Akzeptanzrate"
plot "messen.txt" u 1:3 lt 7 w lines
set xlabel "Temperatur"
set xrange [0:10.5]
plot "messenmittel.txt" u 1:2 lt 7 
unset yrange 
unset xrange
