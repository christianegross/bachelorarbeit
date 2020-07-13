//Christiane, 03.07.20
//Bestimmung der Zeile mit maximalem Durchschnittswert zur Skalierung mit MPI


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "auswertungsfunktionen.h"

int main(int argc, char **argv){
int laengen[3]={12, 24, 48};//Laengen, bei denen gemessen wurde
int prozessorenzahl[2]={1,2};
//benoetigte Variablen
char dateinameskalierung[100];
FILE *dateiskalierung;
int laenge, proz;
//alle prozessorzahlen und laengen
for (int p=0;p<2;p+=1){
	proz=prozessorenzahl[p];
	for (int l=0;l<3;l+=1){
		laenge=laengen[l];
		//dateinamen festlegen
		sprintf(dateinameskalierung, "Messungen/Zeiten/zmittel-m001000-proz00-%.2dmpivmt05-l%d.txt", proz, laenge);
		
		//dateien oeffnen
		dateiskalierung=fopen(dateinameskalierung, "r+");
		maxline(dateiskalierung, 8,4,proz+1);
		fclose(dateiskalierung);
	}
}
return 0;
}
