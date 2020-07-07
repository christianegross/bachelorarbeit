//Christiane, 03.07.20
//Bestimmung einer Zeile mit minimalen Wert fuer Ableitungsminimu/kritischer Punkt


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "auswertungsfunktionen.h"

int main(int argc, char **argv){
int laengen[6]={12, 24, 36, 48, 96, 192};
char dateinameableitung[100], dateinamedreipunkt[100], dateinamebootstrap[100], dateinamemagquad[100], dateinamemagquadsqrt[100];
FILE *dateiableitung, *dateidreipunkt, *dateibootstrap, *dateimagquad, *dateimagquadsqrt;
int laenge;
int temperaturzahl=210, schritt=1;
for (int l=0;l<6;l+=1){
	laenge=laengen[l];
	sprintf(dateinameableitung, "Messungen/ableitung-magnetisierung-laenge-%.4d-m-010240-node02.txt", laenge);
	sprintf(dateinamedreipunkt, "Messungen/ableitungdreipunkt-magnetisierung-laenge-%.4d-m-010240-node02.txt", laenge);
	sprintf(dateinamebootstrap, "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l%.4d-m-010240-node02.txt", laenge);
	sprintf(dateinamemagquad, "Messungen/Bootstrapges/bootstrapalle-magquad-l%.4d-m-010240-node02.txt", laenge);
	sprintf(dateinamemagquadsqrt, "Messungen/Bootstrapges/bootstrapalle-magquadsqrt-l%.4d-m-010240-node02.txt", laenge);
	dateiableitung=fopen(dateinameableitung, "w+");
	dateidreipunkt=fopen(dateinamedreipunkt, "w+");
	dateibootstrap=fopen(dateinamebootstrap, "r");
	dateimagquad=fopen(dateinamemagquad, "r");
	dateimagquadsqrt=fopen(dateinamemagquadsqrt, "w+");
	sqrtspalte(dateimagquad, dateimagquadsqrt, 6,3,10*temperaturzahl);
	ableitung(128, temperaturzahl/schritt*12, 6, 3,4,5,1, dateibootstrap, dateiableitung);
	ableitungdreipunkt(128, temperaturzahl/schritt*12, 6, 3,4,5,1, dateibootstrap, dateidreipunkt);
	//ableitung(128, temperaturzahl/schritt*12, 6, 3,4,5,1, dateimagquadsqrt, dateiableitung);
	//ableitungdreipunkt(128, temperaturzahl/schritt*12, 6, 3,4,5,1, dateimagquadsqrt, dateidreipunkt);
	printf("%e\t", (double)laenge);
	minline(dateiableitung, 3, 1, 209);
	//minline(dateidreipunkt, 3, 1, 208);
	fclose(dateiableitung);
	fclose(dateidreipunkt);
	fclose(dateibootstrap);
	fclose(dateimagquadsqrt);
	fclose(dateimagquad);
}
return 0;
}

