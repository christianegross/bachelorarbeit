//Christiane, 03.07.20
//Bestimmung der Zeile mit maximalem Durchschnittswert zur Skalierung der messzeit bei MPI
//Ausgabe der Zeile, zusaetzlich Bestimmung des speedups


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "auswertungsfunktionen.h"

int main(int argc, char **argv){
//int laengen[5]={12, 60, 240, 840, 1260};//{12, 24, 36, 48, 60, 72, 96, 120, 144, 180, 240, 360, 420, 840, 1260};//Laengen, bei denen gemessen wurde
double temperaturen[4]={0.5, 2.2, 2.3, 10};//Temperaturen, bei denen gemessen wurde
//int prozessorenzahl[5]={1,2,3,4,6};
//benoetigte Variablen
char dateinameskalierung[100];
FILE *dateiskalierung;
int laenge;
double ergebnisse[2], zeiteinproz, vareinproz, temperatur;
//alle prozessorzahlen und laengen

for (int l=0;l<4;l+=1){
	//Auswahl moeglich, ob bei verschiedenen Laengen oder verschiedenen Temperaturen Wert bestimmt wird
	laenge=480;//laengen[l];
	temperatur=temperaturen[l];
	for (int p=1;p<=20;p+=1){//Alle moeglichen Prozesszahlen durchgehen, je Prozesszahl durchschnittliche Laufzeit des langsamsten Prozesses herausfinden und daraus Speedup bestimmen 
		if ((laenge%p==0)&&(laenge!=p)){//Nur machen, wenn Prozessorzahl möglich ist
		//proz=prozessorenzahl[p];
		//dateinamen festlegen
		//sprintf(dateinameskalierung, "Messungen/Zeiten/zmittel-m001000-proz00-0%dmpiskalt05n%d-l%d.txt", proz, proz, laenge);//verschiedene Laengen bei gleicher Temperatur
			sprintf(dateinameskalierung, "Messungen/Zeiten/zmittel-m001000-proz00-%.2dskalmpithermrichtigt%.2d-l%d.txt", p, (int)(10*temperatur), laenge);//verschiedene Temperaturen bei gleicher Laenge
		
			//dateien oeffnen
			dateiskalierung=fopen(dateinameskalierung, "r+");
			maxline(dateiskalierung, 8,4,p+1, ergebnisse, 4, 5);//Maximale Laufzeit bestimmen
			//Aus Ergebnissen Speedup berechnen
			if(p==1){//Zeiten von einem Prozess seperat zum Vergleich benoetigt
				zeiteinproz=ergebnisse[0];
				vareinproz=ergebnisse[1];
				printf("%f\t%f\n", 1.0, vareinproz/zeiteinproz);
			}
			else{
				printf("%f\t%f\n", zeiteinproz/ergebnisse[0], sqrt((vareinproz/ergebnisse[0])*(vareinproz/ergebnisse[0])+(zeiteinproz*ergebnisse[1]/ergebnisse[0]/ergebnisse[0])*(zeiteinproz*ergebnisse[1]/ergebnisse[0]/ergebnisse[0])));
			}
			fclose(dateiskalierung);
		}
	}
	printf("\n");//optisch schoener, Trennung zwischen Bloecken
}
return 0;
}
