//Christiane, 03.07.20
//Bestimmung einer Zeile mit minimalen Wert fuer Ableitungsminimum/kritischer Punkt


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "auswertungsfunktionen.h"

int main(int argc, char **argv){
//~ int laengen[25]={12, 14, 16, 18, 20, 22, 24, 28, 32, 36, 40, 44, 48, 56, 64, 72, 80, 88, 96, 112, 128, 144, 192, 384, 768};//Laengen, bei denen gemessen wurde
int laengen[15]={12, 14, 16, 18, 20, 22, 24, 28, 32, 40, 48, 72, 96, 240, 480};//Laengen, bei denen gemessen wurde
int schritte[1]={1/*, 2, 5, 10*/};//SChritte, mit denen gemessen wurde
//benoetigte Variablen
char dateinameableitung[200], dateinamedreipunkt[200], dateinamebootstrap[200];//, dateinamemagquad[100], dateinamemagquadsqrt[100];
FILE *dateiableitung, *dateidreipunkt, *dateibootstrap, *dateimagquad, *dateimagquadsqrt;
int laenge, schritt, p=2;
int temperaturzahl=210;
//alle schritte und laengen werden durchgegangen
for (int s=0;s<1;s+=1){
	schritt=schritte[s];
	for (int l=0;l<15;l+=1){
		laenge=laengen[l];
		if(laenge>=240){p=20;}
		if((laenge>=72)&&(laenge<240)){p=8;}
		if((laenge>=48)&&(laenge<72)){p=6;}
		if((laenge>=28)&&(laenge<48)){p=4;}
		if((laenge>=10)&&(laenge<28)){p=2;}
		//dateinamen festlegen
		//~ sprintf(dateinameableitung, "Messungen/Ableitungen/ableitung-magnetisierung-laenge-%.4d-m-010240-node02-sch-%.2d.txt", laenge, schritt);
		//~ sprintf(dateinamedreipunkt, "Messungen/Ableitungen/ableitungdreipunkt-magnetisierung-laenge-%.4d-m-010240-node02-sch-%.2d.txt", laenge, schritt);
		//~ sprintf(dateinamebootstrap, "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l%.4d-m-010240-node00-sch-%.2d.txt", laenge, schritt);
		sprintf(dateinameableitung, "Messungen/Ableitungen/ableitung-magnetisierung-laenge-%.4d-m-010000-node02-proz%.2d-sch-%.2d.txt", laenge, p, schritt);
		sprintf(dateinamedreipunkt, "Messungen/Ableitungen/ableitungdreipunkt-magnetisierung-laenge-%.4d-m-010000-node02-proz%.2d-sch-%.2d.txt", laenge, p, schritt);
		sprintf(dateinamebootstrap, "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l%.4d-m-010000-proz%.2d-sch-%.2d.txt", laenge, p, schritt);
		//sprintf(dateinamebootstrap, "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l0022-m-010000-proz02-sch-01.txt", laenge, p, schritt);
		//dateien oeffnen
		//dateien, in die die Ableitungen gespeichert werden sollen
		dateiableitung=fopen(dateinameableitung, "w+");//zweipunktmethode
		dateidreipunkt=fopen(dateinamedreipunkt, "w+");//dreipunktmethode
		dateibootstrap=fopen(dateinamebootstrap, "r");//datei mit den Messwerten des Betrags der Magnetisierung
		//bestimmt Ableitungen nach zwei und dreipunktmethode der einfachen Magnetisierungswerte
		ableitung(128, temperaturzahl/schritt*12, 6, 3,4,5,1, dateibootstrap, dateiableitung);
		ableitungdreipunkt(128, temperaturzahl/schritt*12, 6, 3,4,5,1, dateibootstrap, dateidreipunkt);
		//schreibt ob zwei oder dreipunkt, schritt, laenge in Ausgabe, Inhalt der minimalen zeile aus minline, 2, 3 fuer Ableitung, 0,1 fuer norm/sqrt
		printf("%f\t%f\t2.0\t0.0\t", (double)schritt, (double)laenge);
		minline(dateiableitung, 3, 1, 209); //bestimmt minimale Zeile aus zweipunktableitung und gibt sie aus
		printf("%f\t%f\t3.0\t0.0\t", (double)schritt, (double)laenge);
		minline(dateidreipunkt, 3, 1, 208);//bestimmt minimale Zeile aus dreipunktableitung und gibt sie aus
		//close, damit keine Speicherprobleme entstehen
		fclose(dateiableitung);
		fclose(dateidreipunkt);
		fclose(dateibootstrap);
	}
}
return 0;
}

