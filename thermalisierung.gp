set ter pdfcairo size 8in, 11in
set out "thermalisierungplot.pdf"

plot "thermalisierung.txt" u 1:2 lt 7 w lines

set yrange [0:200*200]
plot "messen.txt" u 1:2 lt 7 w lines
unset yrange 
