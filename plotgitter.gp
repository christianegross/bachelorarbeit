set ter pdfcairo size 11in,8in

set view map
#set dgrid3d


set out "plotgitter.pdf"
set title "Thermalisiert" font "arial,40"
#set cbrange[-1:1]
splot "gitterthermalisiert.txt" using 1:2:3 w image
