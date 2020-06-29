#Variablen
LIBS:= -lgsl -lgslcblas -lm -fopenmp -lrt

#auswertung: auswertungzeit.o
#	gcc -std=c99 -Wall -pedantic -o $@ $^ $(LIBS)

#skalfak00: skalierungfak.o
#	gcc -std=c99 -Wall -pedantic -o $@ $^ $(LIBS)

skalierungnode00: skalierung.o messfunktionen.o sweeps.o auswertungsfunktionen.o 
	gcc -std=c99 -Wall -pedantic -o $@ $^ $(LIBS)

#skalierung1606: skalierung.o messfunktionen1606.o auswertungsfunktionen.o sweeps1606.o 
#	gcc -std=c99 -Wall -pedantic -o $@ $^ $(LIBS)
#sweeps1606.o 

#erstellt programm ising
#ising: ising.o messfunktionen.o sweeps.o auswertungsfunktionen.o
#	gcc -std=c99 -Wall -pedantic -o $@ $^ $(LIBS)


#bootstrap: bootstrap.o messfunktionen.o sweeps.o auswertungsfunktionen.o
#	gcc -std=c99 -Wall -pedantic -o $@ $^ $(LIBS)
	
#erstellt aus allen .c dateien eine .o datei	
%.o: %.c
	gcc -std=c99 -Wall -pedantic -fopenmp -I /usr/include/ $^ -c 
	
#.PHONY: clean

#clean:

#$^: fügt alle abhängigkeiten ein
#%: sucht alle Dateien mit gegebener Endung(wildcard)
#$@: Fügt target ein
