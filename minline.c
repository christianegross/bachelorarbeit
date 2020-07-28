//Christiane, 03.07.20
//Bestimmung einer Zeile mit minimalen Wert fuer Ableitungsminimum/kritischer Punkt


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "auswertungsfunktionen.h"

int main(int argc, char **argv){
int laengen[25]={12, 14, 16, 18, 20, 22, 24, 28, 32, 36, 40, 44, 48, 56, 64, 72, 80, 88, 96, 112, 128, 144, 192, 384, 768};//Laengen, bei denen gemessen wurde
int schritte[1]={1/*, 2, 5, 10*/};//SChritte, mit denen gemessen wurde
//benoetigte Variablen
char dateinameableitung[100], dateinamedreipunkt[100], dateinamebootstrap[100], dateinamemagquad[100], dateinamemagquadsqrt[100];
FILE *dateiableitung, *dateidreipunkt, *dateibootstrap, *dateimagquad, *dateimagquadsqrt;
int laenge, schritt;
int temperaturzahl=210;
//alle schritte und laengen werden durchgegangen
for (int s=0;s<1;s+=1){
	schritt=schritte[s];
	for (int l=0;l<25;l+=1){
		laenge=laengen[l];
		//dateinamen festlegen
		sprintf(dateinameableitung, "Messungen/Ableitungen/ableitung-magnetisierung-laenge-%.4d-m-010368-node02-sch-%.2d.txt", laenge, schritt);
		sprintf(dateinamedreipunkt, "Messungen/Ableitungen/ableitungdreipunkt-magnetisierung-laenge-%.4d-m-010368-node02-sch-%.2d.txt", laenge, schritt);
		sprintf(dateinamebootstrap, "Messungen/Bootstrapges/bootstrapalle-magnetisierung-l%.4d-m-010368-node02-sch-%.2d.txt", laenge, schritt);
		sprintf(dateinamemagquad, "Messungen/Bootstrapges/bootstrapalle-magquad-l%.4d-m-010368-node02-sch-%.2d.txt", laenge, schritt);
		sprintf(dateinamemagquadsqrt, "Messungen/Bootstrapges/bootstrapalle-magquadsqrt-l%.4d-m-010368-node02-sch-%.2d.txt", laenge, schritt);
		//dateien oeffnen
		//dateien, in die die Ableitungen gespeichert werden sollen
		dateiableitung=fopen(dateinameableitung, "w+");//zweipunktmethode
		dateidreipunkt=fopen(dateinamedreipunkt, "w+");//dreipunktmethode
		dateibootstrap=fopen(dateinamebootstrap, "r");//datei mit den Messwerten des Betrags der Magnetisierung
		dateimagquad=fopen(dateinamemagquad, "r");//datei mit den Messwerten des Quadrats der Magnetisierung
		dateimagquadsqrt=fopen(dateinamemagquadsqrt, "w+");//datei, in die Wurzel aus Quadrat der Magnetisierung gespeichert wird
		sqrtspalte(dateimagquad, dateimagquadsqrt, 6,3,10*temperaturzahl);//bestimmt Wurzel aus Quadrat der Magnetisierung, gibt alle anderen Spalten unveraendert aus
		//bestimmt Ableitungen nach zwei und dreipunktmethode der einfachen Magnetisierungswerte
		ableitung(256, temperaturzahl/schritt*12, 6, 3,4,5,1, dateibootstrap, dateiableitung);
		ableitungdreipunkt(256, temperaturzahl/schritt*12, 6, 3,4,5,1, dateibootstrap, dateidreipunkt);
		//schreibt ob zwei oder dreipunkt, schritt, laenge in Ausgabe, Inhalt der minimalen zeile aus minline, 2, 3 fuer Ableitung, 0,1 fuer norm/sqrt
		printf("%f\t%f\t2.0\t0.0\t", (double)schritt, (double)laenge);
		minline(dateiableitung, 3, 1, 209); //bestimmt minimale Zeile aus zweipunktableitung und gibt sie aus
		printf("%f\t%f\t3.0\t0.0\t", (double)schritt, (double)laenge);
		minline(dateidreipunkt, 3, 1, 208);//bestimmt minimale Zeile aus dreipunktableitung und gibt sie aus
		//Schliesst und oeffnet->ueberschreibt Ableitungsdateien, damit minline korrekt die zweiten bestimmten Ableitungen betrachtet
		//Ableitungen nicht wichtig, nur minimum
		fclose(dateiableitung);
		fclose(dateidreipunkt);
		dateiableitung=fopen(dateinameableitung, "w+");//zweipunktmethode
		dateidreipunkt=fopen(dateinamedreipunkt, "w+");//dreipunktmethode
		//bestimmt Ableitungen nach zwei und dreipunktmethode der Wurzel aus dem Quadrat der Magnetisierung
		ableitung(256, temperaturzahl/schritt*12, 6, 3,4,5,1, dateimagquadsqrt, dateiableitung);
		ableitungdreipunkt(256, temperaturzahl/schritt*12, 6, 3,4,5,1, dateimagquadsqrt, dateidreipunkt);
		//schreibt ob zwei oder dreipunkt, schritt, laenge in Ausgabe, Inhalt der minimalen zeile aus minline
		printf("%f\t%f\t2.0\t1.0\t", (double)schritt, (double)laenge);
		minline(dateiableitung, 3, 1, 209); //bestimmt minimale Zeile aus zweipunktableitung und gibt sie aus
		printf("%f\t%f\t3.0\t1.0\t", (double)schritt, (double)laenge);
		minline(dateidreipunkt, 3, 1, 208);//bestimmt minimale Zeile aus dreipunktableitung und gibt sie aus
		//close, damit keine Speicherprobleme entstehen
		fclose(dateiableitung);
		fclose(dateidreipunkt);
		fclose(dateibootstrap);
		fclose(dateimagquadsqrt);
		fclose(dateimagquad);
	}
}
return 0;
}

