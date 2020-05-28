//Christiane, 25.05.20
//Skalierungstest von sweep fuer verschiedene cores


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "messfunktionen.h"
#include "auswertungsfunktionen.h"

int main(int argc, char **argv){
	//benoetigte Variablen initialisieren
	int maxcores=omp_get_max_threads();//aus Computerarchitektut/batchskript
	int laenge=10;//laenge der verwendeten Gitter
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int messungen=1000;//pro temperatur
	double temperatur=3.5;//Skalierung bei nur einer TEmperatur messen
	FILE *messdatei, *dummydatei, *zeitdatei;//speichern der Messergenbnisse der einzelnen Messungen, der Zeiten und ein nicht benoetigtes nicht thermalisiertes Gitter
	char dateinamemessen[150], dateinamezeit[150];
	struct timeval anfangmessen, endemessen;//Zeitmessung mit gettimeofday
	double sec, usec, zeitgesmessen, zeiteincore;
	int gitter[laenge*laenge];//Gitter erstellen und von thermalisieren ausgeben lassen
	initialisierung(gitter, laenge, seed);
	dummydatei=fopen("dummytherm.txt", "w+");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird	
	sprintf(dateinamezeit,"Messungen/Zeiten/zeitenmessen-laenge%.4d-m%.6d-mehrere.txt",laenge,messungen);
	zeitdatei=fopen(dateinamezeit, "w+");
	gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	gsl_rng_set(generator, seed);
	thermalisieren(laenge, temperatur, j, seed, 1, gitter, dummydatei, generator);//Erstes Thermalisieren, hier nur zur Ausgabe des Gitters
	for (int durchlauf=0; durchlauf<100;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden
	//Vergleichsmassstab: Messungen bei einem core	
		gsl_rng_set(generator, seed);
		double speedup=1;//aus definition
		omp_set_num_threads(1);
		sprintf(dateinamemessen,"Messungen/Messwerte/messung-laenge%.4d-m%.6d-cores%.2d-%.2d.txt",laenge,messungen,1, durchlauf);
		messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
		gettimeofday(&anfangmessen, NULL);
		messen(laenge, temperatur, j, messungen, dummydatei, messdatei, generator);
		gettimeofday(&endemessen, NULL);
		sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
		usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
		zeiteincore=sec+1e-06*usec;
		printf("bei %d cores haben %d Messungen %f Sekunden gebraucht\n", 1, messungen, zeiteincore);
		fprintf(zeitdatei, "%f\t%f\t%f\t%f\n", 1.0, (double)messungen, zeiteincore, 1.0);//cores messungen Zeit Speedup
		fclose(messdatei);
		//Messungen bei mehreren cores
		for (int cores=2;cores<=maxcores;cores+=1){
			omp_set_num_threads(cores);
			sprintf(dateinamemessen,"Messungen/Messwerte/messung-laenge%.4d-m%.6d-cores%.2d-%.2d.txt",laenge,messungen,cores,durchlauf);
			messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
			gettimeofday(&anfangmessen, NULL);
			messen(laenge, temperatur, j, messungen, dummydatei, messdatei, generator);
			gettimeofday(&endemessen, NULL);
			sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
			usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
			zeitgesmessen=sec+1e-06*usec;
			speedup=zeitgesmessen/zeiteincore;
			printf("bei %d cores haben %d Messungen %f Sekunden gebraucht\n", cores, messungen, zeitgesmessen);
			fprintf(zeitdatei, "%f\t%f\t%f\t%f\n", (double)cores, (double)messungen, zeitgesmessen, speedup);
			fclose(messdatei);
		}
	}
	fclose(dummydatei);
	
	gsl_rng_free(generator);//free, close: zum Verhindern von Speicherproblemen
}
