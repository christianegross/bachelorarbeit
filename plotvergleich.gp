set ter pdfcairo size 8in,5.5in


#Magnetisierung
set out "Messungen/Vergleich.pdf"

set style line 1 lt 7 lc 1
set style line 2 lt 7 lc 2
set style line 3 lt 7 lc 3
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

do for [t=0:299:30]{
datei1=sprintf("Messungen/Vergleichwerte/vergleich-laenge0050-t%.3d.txt", t)
datei2=sprintf("Messungen/Messwerte/messung-laenge0050-t%.3d.txt", t)

set title sprintf("T=%f, Hamiltonian", t*0.015+0.015)
plot datei1 u 1:3 ls 3 ps 0.2 title title3,\
	datei1 u 1:4 ls 4 ps 0.2 title title4,\
	datei1 u 1:7 ls 7 ps 0.2 title title7

set title sprintf("T=%f, Magnetisierung", t*0.015+0.015)
#plot datei2 u 1:5 ls 3 ps 0.2 title title3,\
#	datei2 u 1:7 ls 4 ps 0.2 title title4,\
#	datei2 u 1:($5-$7) ls 7 ps 0.2 title title7,\

}
