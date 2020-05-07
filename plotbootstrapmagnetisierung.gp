set ter pdfcairo size 4in,5.5in

set out "Messungen/magnetisierungmitbootstrap.pdf"

do for [n in "8 16 32 64 128 256 384 512 640 758 876 1024"]{
set title "".n
plot "Messungen/bootstrapalle-l0050-m-016384.txt" using (($1==n)?$5:1/0):3:4 w yerrorbars lt 7 ps 0.3 title ""

}

do for [n in "32 64 128 256 384 512 640 758 876 1024 1280 1536"]{
set title "".n
#plot "Messungen/bootstrapalle-l0050-m-032768.txt" using (($1==n)?$5:1/0):3:4 w yerrorbars lt 7 ps 0.3 title ""

}
