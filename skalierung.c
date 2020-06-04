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
	int lenarray[8]={6, 10,20, 36, 50, 100, 200, 350};
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int messungen=1000;//pro temperatur
	double mittelzeit, varianzzeit, speedupmittel, speedupfehler, speedup;
	int durchlaeufe=10;
	double *ergebnisse;
	if((ergebnisse=(double*)malloc(sizeof(double)*durchlaeufe))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren der Temperaturen!\n");
		return (-1);
	}
	//double mittelham, varianzham;
	double temperatur=3.5;//Skalierung bei nur einer TEmperatur messen
	FILE *messdatei, *dummydatei, *zeitdatei, *mitteldatei;//speichern der Messergenbnisse der einzelnen Messungen, der Zeiten und ein nicht benoetigtes nicht thermalisiertes Gitter
	char dateinamezeit[200], dateinamemittel[200];//dateinamemessen[150],
	struct timeval anfangmessen, endemessen;//Zeitmessung mit gettimeofday
	double sec, usec, zeitgesmessen, zeiteincore, varianzeincore;
	dummydatei=fopen("dummytherm.txt", "w+");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird	
	//sprintf(dateinamezeit,"Messungen/Zeiten/zeitenmessen-laenge%.4d-m%.6d-mehrere.txt",laenge,messungen);
	sprintf(dateinamezeit,"Messungen/Zeiten/zeitenmessen-m%.6d-mehrerelaengenunddurchlaeufenode01dynamic.txt",messungen);
	sprintf(dateinamemittel,"Messungen/Zeiten/zeitenmittel-m%.6d-mehrerelaengenunddurchlaeufenode01dynamic.txt",messungen);
	zeitdatei=fopen(dateinamezeit, "w+");
	mitteldatei=fopen(dateinamemittel, "w+");
	gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	gsl_rng_set(generator, seed);
	for (int laengen=0; laengen<8; laengen+=1){
		laenge=lenarray[laengen];
		printf("Laenge=%d\n", laenge);
		int gitter[laenge*laenge];//Gitter erstellen und von thermalisieren ausgeben lassen
		initialisierung(gitter, laenge, seed);
		thermalisieren(laenge, temperatur, j, seed, 1, gitter, dummydatei, generator);//Erstes Thermalisieren, hier nur zur Ausgabe des Gitters
		for (int durchlauf=0; durchlauf<durchlaeufe;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden
		//Vergleichsmassstab: Messungen bei einem core	
			gsl_rng_set(generator, seed);
			speedup=1;//aus definition
			omp_set_num_threads(1);
			//sprintf(dateinamemessen,"Messungen/Messwerte/messung-laenge%.4d-m%.6d-cores%.2d-%.2d.txt",laenge,messungen,1, durchlauf);
			//messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
			messdatei=fopen("dummymessen.txt", "w+");
			//vergleichsdatei=fopen("dummyvergleich.txt", "w+");
			gettimeofday(&anfangmessen, NULL);
			//messenvergleichen(laenge, temperatur, j, messungen, dummydatei, messdatei, vergleichsdatei, generator);
			messen(laenge, temperatur, j, messungen, dummydatei, messdatei, generator);
			gettimeofday(&endemessen, NULL);
			//~ mittelham=mittelwertberechnungnaiv(vergleichsdatei, messungen, 3, 4);
			//~ varianzham=varianzberechnungnaiv(vergleichsdatei, messungen, mittelham, 3, 4);
			//~ printf("Laenge %d, H=%f+/-%f\n", laenge, mittelham, varianzham);
			sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
			usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
			zeiteincore=sec+1e-06*usec;
			ergebnisse[durchlauf]=zeiteincore;
			printf("bei %d cores haben %d Messungen %f Sekunden gebraucht\n", 1, messungen, zeiteincore);
			fprintf(zeitdatei, "%f\t%f\t%f\t%f\t%f\n", 1.0, (double)messungen, zeiteincore, 1.0, (double)laenge);//cores messungen Zeit Speedup
			fclose(messdatei);
		}
		mittelzeit=mittelwertarray(ergebnisse, durchlaeufe);
		varianzzeit=varianzarray(ergebnisse, durchlaeufe, mittelzeit);
		zeiteincore=mittelzeit;
		varianzeincore=varianzzeit;
		speedupmittel=1;
		speedupfehler=varianzzeit/mittelzeit;
		fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\n", 1.0, (double)laenge, mittelzeit, varianzzeit, speedupmittel, speedupfehler);
			//fclose(vergleichsdatei);
			//Messungen bei mehreren cores
		for (int cores=2;cores<=maxcores;cores+=1){
			for (int durchlauf=0; durchlauf<durchlaeufe;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden
				omp_set_num_threads(cores);
				//sprintf(dateinamemessen,"Messungen/Messwerte/messung-laenge%.4d-m%.6d-cores%.2d-%.2d.txt",laenge,messungen,cores,durchlauf);
				//messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
				messdatei=fopen("dummymessen.txt", "w+");
				//vergleichsdatei=fopen("dummyvergleich.txt", "w+");
				gettimeofday(&anfangmessen, NULL);
				//messenvergleichen(laenge, temperatur, j, messungen, dummydatei, messdatei, vergleichsdatei, generator);
				messen(laenge, temperatur, j, messungen, dummydatei, messdatei, generator);
				gettimeofday(&endemessen, NULL);
				//~ mittelham=mittelwertberechnungnaiv(vergleichsdatei, messungen, 3, 4);
				//~ varianzham=varianzberechnungnaiv(vergleichsdatei, messungen, mittelham, 3, 4);
				//~ printf("Laenge %d, H=%f+/-%f\n", laenge, mittelham, varianzham);
				sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
				usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
				zeitgesmessen=sec+1e-06*usec;
				ergebnisse[durchlauf]=zeitgesmessen;
				speedup=zeitgesmessen/zeiteincore;
				printf("bei %d cores haben %d Messungen %f Sekunden gebraucht\n", cores, messungen, zeitgesmessen);
				fprintf(zeitdatei, "%f\t%f\t%f\t%f\t%f\n", (double)cores, (double)messungen, zeitgesmessen, speedup, (double)laenge);
				fclose(messdatei);
		//fclose(vergleichsdatei);
			}
			mittelzeit=mittelwertarray(ergebnisse, durchlaeufe);
			varianzzeit=varianzarray(ergebnisse, durchlaeufe, mittelzeit);
			speedupmittel=mittelzeit/zeiteincore;
			speedupfehler=sqrt(((varianzzeit/zeiteincore)*(varianzzeit/zeiteincore))+(mittelzeit*varianzeincore/zeiteincore/zeiteincore)*(mittelzeit*varianzeincore/zeiteincore/zeiteincore));
			fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\n", (double)cores, (double)laenge, mittelzeit, varianzzeit, speedupmittel, speedupfehler);
		}
	}
	fclose(dummydatei);
	fclose(zeitdatei);
	
	gsl_rng_free(generator);//free, close: zum Verhindern von Speicherproblemen
	
	//free(ergebnisse);
}
