#plottet hamiltonian und akzeptanzrate, die durch unterschiedliche sweep-Funktionen gemessen wurde, plottet teilweide auch deren Differenz
set terminal pdfcairo size 4in, 2.75in color colortext


#Magnetisierung
set output "Vergleich.pdf"

set style line 1 lc 1 pt 4
set style line 2 lc 2 pt 8 
set style line 3 lc 3 pt 10
set style line 4 lt 7 lc 4
set style line 5 lt 7 lc 5
set style line 6 lt 7 lc 6
set style line 7 lt 7 lc 7
set style line 8 lt 7 lc 8


title5="alt-ohnepar"
title6="alt-mitpar"
title7="ohnepar-mitpar"
title2="alt"
title3="ohnepar"
title4="mitpar"

set title "Hamiltonian"
set xlabel "Temperatur"
set ylabel "H"
set key bottom right font 'arial, 12'

dateimittel="Messungen/Mittelwerte/vergleichmittel-l0050-m-010000.txt"

plot dateimittel u 1:2:3 w yerrorbars ls 1 ps 0.4 title "Zeilenweise",\
dateimittel u 1:4:5 w yerrorbars ls 2 ps 0.4 title "Schachbrett nicht parallel",\
dateimittel u 1:6:7 w yerrorbars ls 3 ps 0.4 title "Schachbrett parallel",\

set title "Akzeptanzrate"
set ylabel "Akzeptanzrate"

plot dateimittel u 1:8:9 w yerrorbars ls 1 ps 0.4 title "Zeilenweise",\
dateimittel u 1:10:11 w yerrorbars ls 2 ps 0.4 title "Schachbrett nicht parallel",\
dateimittel u 1:12:13 w yerrorbars ls 3 ps 0.4 title "Schachbrett parallel",\

#do for [t=0:299:30]{
#datei1=sprintf("Messungen/Vergleichwerte/vergleich-laenge0050-t%.3d.txt", t)
#datei2=sprintf("Messungen/Messwerte/messung-laenge0050-t%.3d.txt", t)

#set title sprintf("T=%f, Hamiltonian", t*0.015+0.015)
#plot datei1 u 1:3 ls 3 ps 0.4 title title3,\
#	datei1 u 1:4 ls 4 ps 0.4 title title4,\
#	datei1 u 1:7 ls 7 ps 0.4 title title7
#}
#do for [t=0:299:30]{
#datei1=sprintf("Messungen/Vergleichwerte/vergleich-laenge0050-t%.3d.txt", t)
#datei2=sprintf("Messungen/Messwerte/messung-laenge0050-t%.3d.txt", t)

#set title sprintf("T=%f, Magnetisierung", t*0.015+0.015)
#plot datei2 u 1:5 ls 3 ps 0.4 title title3,\
#	datei2 u 1:7 ls 4 ps 0.4 title title4,\
#	datei2 u 1:($5-$7) ls 7 ps 0.4 title title7,\

#}
