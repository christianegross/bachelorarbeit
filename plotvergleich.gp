set ter pdfcairo size 4in,5.5in


#Magnetisierung
set out "Messungen/Vergleich.pdf"

#plot "Messungen/Vergleichwerte/vergleich-laenge0050-t250.txt" u 1:2 lt 7 ps 0.3 title "Unterschied bei verschiedenen sweep-arten"

set title "T 250, Hamiltonian"
plot "Messungen/Vergleichwerte/vergleich-laenge0050-t250.txt" u 1:2 ps 0.1 lt 7 title "Unterschied bei verschiedenen sweep-arten", "Messungen/Vergleichwerte/vergleich-laenge0050-t250.txt" u 1:3 ps 0.1 lt 1 title "H", "Messungen/Vergleichwerte/vergleich-laenge0050-t250.txt" u 1:4 ps 0.1 lt 6 title "Halt"

set title "T 250, Akzeptanzrate"
plot "Messungen/Messwerte/messung-laenge0050-t250.txt" u 1:2 ps 0.1 lt 1 title "neu", "Messungen/Messwerte/messung-laenge0050-t250.txt" u 1:4 ps 0.1 lt 6 title "alt","Messungen/Messwerte/messung-laenge0050-t250.txt" u 1:($2-$4) ps 0.1 lt 7 title "differenz"

set title "T 150, Hamiltonian"
plot "Messungen/Vergleichwerte/vergleich-laenge0050-t150.txt" u 1:2 ps 0.1 lt 7 title "Unterschied bei verschiedenen sweep-arten", "Messungen/Vergleichwerte/vergleich-laenge0050-t150.txt" u 1:3 ps 0.1 lt 1 title "H", "Messungen/Vergleichwerte/vergleich-laenge0050-t150.txt" u 1:4 ps 0.1 lt 6 title "Halt"

set title "T 150, Akzeptanzrate"
plot "Messungen/Messwerte/messung-laenge0050-t150.txt" u 1:2 ps 0.1 lt 1 title "neu", "Messungen/Messwerte/messung-laenge0050-t150.txt" u 1:4 ps 0.1 lt 6 title "alt","Messungen/Messwerte/messung-laenge0050-t150.txt" u 1:($2-$4) ps 0.1 lt 7 title "differenz"

set title "T 50, Hamiltonian"
plot "Messungen/Vergleichwerte/vergleich-laenge0050-t050.txt" u 1:2 ps 0.1 lt 7 title "Unterschied bei verschiedenen sweep-arten", "Messungen/Vergleichwerte/vergleich-laenge0050-t050.txt" u 1:3 ps 0.1 lt 1 title "H", "Messungen/Vergleichwerte/vergleich-laenge0050-t050.txt" u 1:4 ps 0.1 lt 6 title "Halt"

set title "T 50, Akzeptanzrate"
#set yrange[-0.01:0.01]
plot "Messungen/Messwerte/messung-laenge0050-t050.txt" u 1:2 ps 0.1 lt 1 title "neu", "Messungen/Messwerte/messung-laenge0050-t050.txt" u 1:4 ps 0.1 lt 6 title "alt","Messungen/Messwerte/messung-laenge0050-t050.txt" u 1:($2-$4) ps 0.1 lt 7 title "differenz"
