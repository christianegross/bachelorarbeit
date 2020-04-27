set ter pdfcairo size 8in,8in

set view map
#set hidden3d


set out "Messungen/plotgitter.pdf"
set palette defined (0 "white" , 50 "gray", 100 "black")
set title "Thermalisiert" font "arial,40"
set cbrange[-1:1]

#laenge 30
set xrange[-0.5:29.5]
set yrange[-0.5:29.5]

set title "T=0.25, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t00.txt" using 1:2:3 w image title ""

set title "T=0.50, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t01.txt" using 1:2:3 w image title ""

set title "T=0.75, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t02.txt" using 1:2:3 w image title ""

set title "T=1.00, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t03.txt" using 1:2:3 w image title ""

set title "T=1.50, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t04.txt" using 1:2:3 w image title ""

set title "T=2.00, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t05.txt" using 1:2:3 w image title ""

set title "T=2.50, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t06.txt" using 1:2:3 w image title ""

set title "T=3.00, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t07.txt" using 1:2:3 w image title ""

set title "T=4.00, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t08.txt" using 1:2:3 w image title ""

set title "T=5.00, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t09.txt" using 1:2:3 w image title ""

set title "T=6.00, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t10.txt" using 1:2:3 w image title ""

set title "T=7.00, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t11.txt" using 1:2:3 w image title ""

set title "T=8.00, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t12.txt" using 1:2:3 w image title ""

set title "T=9.00, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t13.txt" using 1:2:3 w image title ""

set title "T=10.00, Laenge=0030" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0030-t14.txt" using 1:2:3 w image title ""

unset xrange
unset yrange

#laenge 300
set xrange[-0.5:299.5]
set yrange[-0.5:299.5]

set title "T=0.25, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t00.txt" using 1:2:3 w image title ""

set title "T=0.50, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t01.txt" using 1:2:3 w image title ""

set title "T=0.75, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t02.txt" using 1:2:3 w image title ""

set title "T=1.00, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t03.txt" using 1:2:3 w image title ""

set title "T=1.50, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t04.txt" using 1:2:3 w image title ""

set title "T=2.00, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t05.txt" using 1:2:3 w image title ""

set title "T=2.50, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t06.txt" using 1:2:3 w image title ""

set title "T=3.00, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t07.txt" using 1:2:3 w image title ""

set title "T=4.00, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t08.txt" using 1:2:3 w image title ""

set title "T=5.00, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t09.txt" using 1:2:3 w image title ""

set title "T=6.00, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t10.txt" using 1:2:3 w image title ""

set title "T=7.00, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t11.txt" using 1:2:3 w image title ""

set title "T=8.00, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t12.txt" using 1:2:3 w image title ""

set title "T=9.00, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t13.txt" using 1:2:3 w image title ""

set title "T=10.00, Laenge=0300" font "arial,40"
splot "Messungen/ThermalisierteGitter/thermalisierung-l0300-t14.txt" using 1:2:3 w image title ""

