set ter pdfcairo size 8in,8in

set view map
#set dgrid3d


set out "plotgitter.pdf"
set palette defined (0 "red" , 50 "white", 100 "blue")
set title "Thermalisiert" font "arial,40"
#set cbrange[-1:1]
splot "gitterthermalisiert.txt" using 1:2:3 w image
