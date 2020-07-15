#plottet fuer verschiedene Laenge, Messzahlen Magnetisierung oder ableizung Magnetisierung in Abhängigkeit von der Temperatur
set ter pdfcairo size 4in,5.5in

set out "Messungen/magnetisierungmitbootstrap.pdf"

set xlabel "Temperatur"
set ylabel "Magnetisierung"

do for [n in "8 16 32 64 128 256 384 512 640 758 876 1024"]{
set title "".n
#plot "Messungen/Bootstrapges/bootstrapalle-l0050-m-016384.txt" using (($1==n)?$5:1/0):3:4 w yerrorbars lt 7 ps 0.3 title ""

}

do for [n in "32 64 128 256 384 512 640 758 876 1024 1280 1536"]{
set title "".n
#plot "Messungen/Bootstrapges/bootstrapalle-l0050-m-032768.txt" using (($1==n)?$5:1/0):3:4 w yerrorbars lt 7 ps 0.3 title ""
}

do for [n in "32 64 128 256 384 512 640 758 876 1024 1280 1536"]{
set title "".n
#plot "Messungen/Bootstrapges/bootstrapalle-l0050-m-004096.txt" using (($1==n)?$5:1/0):3:4 w yerrorbars lt 7 ps 0.3 title ""
}


do for [n in "32 64 128 256 384 512 640 758 876 1024 1280 1536"]{
set title "l=".n
plot "Messungen/Bootstrapges/bootstrapalle-l0050-m-010000.txt" using ((($2==n)&&($1==1))?$6:1/0):4:5 w yerrorbars lt 7 ps 0.3 title ""
set title "Ableitung"
x0=NaN
y0=NaN
plot "Messungen/Bootstrapges/bootstrapalle-l0050-m-010000.txt" using ((($2==n)&&($1==1))?(dx=$6-x0,x0=$6,$6-dx/2):1/0):(dy=$4-y0,y0=$4,dy/dx) lt 7 title ""

}

do for [n in "32 64 128 256 384 512 640 758 876 1024 1280 1536"]{
set title "l=".n
plot "Messungen/Bootstrapges/bootstrapalle-l0100-m-010000.txt" using ((($2==n)&&($1==1))?$6:1/0):4:5 w yerrorbars lt 7 ps 0.3 title ""
}

do for [n in "32 64 128 256 384 512 640 758 876 1024 1280 1536"]{
set title "l=".n
plot "Messungen/Bootstrapges/bootstrapalle-l0200-m-010000.txt" using ((($2==n)&&($1==1))?$6:1/0):4:5 w yerrorbars lt 7 ps 0.3 title ""
}

set out "vergleichbootstrap.pdf"
do for [n in "32 64 128 256 384 512 640 758 876 1024 1280 1536"]{
set title "l=".n
#plot "Messungen/Bootstrapges/bootstrapalle-l0050-m-010000.txt" using ((($2==n)&&($1==1))?$6:1/0):4:5 w yerrorbars lt 7 ps 0.3 title "mit par",\
#	"Messungen/Bootstrapges/bootstrapalle-l0050-m-010000.txt" using ((($2==n)&&($1==2))?$6:1/0):4:5 w yerrorbars pt 10 lc 6 ps 0.3 title "ohne par"
}

