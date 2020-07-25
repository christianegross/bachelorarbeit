//Christiane, 25.05.20
//Skalierungstest von sweepmehreregenratoren mit OpenMP fuer verschiedene cores


#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "messfunktionen.h"
#include "sweeps.h"
#include "auswertungsfunktionen.h"


int main(int argc, char **argv){
	//benoetigte Variablen initialisieren
	int maxcores=omp_get_max_threads();//aus Computerarchitektut/batchskript
	int laenge;
	double temperatur;
	char merkmal[50];
	if (argc<3){
		printf("Nicht genug Argumente!\n");
		fprintf(stderr,"Nicht genug Argumente!\n");
		laenge=12;
		temperatur=0.5;
		sprintf(merkmal,"nichtgenugargumente");
	}
	else{	//Parameter koennen aus kommandozeile zugewiesen werden
	laenge=atoi(argv[1]);//laenge der verwendeten Gitter
	temperatur=atof(argv[2]);//Temperatur
	sprintf(merkmal,"%s-l%s",argv[3], argv[1]);
	}
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int messungen=1000;//pro temperatur
	//Variablen, die zum messen benoetigt werden
	double mittelzeit, varianzzeit, speedupmittel, speedupfehler, speedup;//Durschnitt
	double zeitmin, zeitmineincore, speedupmin;//Minimum
	int node=0;//1,2 qbig, 0 vm zur Kennzeichnung
	int durchlaeufe=10;//Anzahl Messungen, ueber die der Durchschnitt gebildet wird
	double *ergebnisse;//Speichert Zeiten der einzelnen Durchlaeufe
	if((ergebnisse=(double*)malloc(sizeof(double)*durchlaeufe))==NULL){//speichert verwendete Temperaturen, 
		fprintf(stderr, "Fehler beim Allokieren der Temperaturen!\n");
		return (-1);
	}
	FILE *messdatei, *dummydatei, *zeitdatei, *mitteldatei;//speichern der Messergenbnisse der einzelnen Messungen, der Zeiten und ein nicht benoetigtesGitter
	char dateinamezeit[200], dateinamemittel[200], dateinamedummytherm[50], dateinamedummymessen[50];//dateinamemessen[150],
	struct timeval anfangmessen, endemessen;//Zeitmessung mit gettimeofday
	double sec, usec, zeitgesmessen, zeiteincore, varianzeincore;
	//Speichert Gitter/Messwerte, diese werden aber nicht benoetigt und daher auch wieder ueberschrieben
	sprintf(dateinamedummytherm, "dummytherm%.2d.txt", node);
	sprintf(dateinamedummymessen, "dummymessen%.2d.txt", node);

	dummydatei=fopen(dateinamedummytherm, "w+");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird	
	sprintf(dateinamezeit,"Messungen/Zeiten/zmessen-m%.6d-node%.2d%s.txt",messungen, node, merkmal);
	sprintf(dateinamemittel,"Messungen/Zeiten/zmittel-m%.6d-node%.2d%s.txt",messungen, node, merkmal);
	zeitdatei=fopen(dateinamezeit, "w+");
	mitteldatei=fopen(dateinamemittel, "w+");
	gsl_rng **generatoren;
	if ((generatoren=(gsl_rng**)malloc(maxcores*sizeof(gsl_rng**)))==NULL){
		printf("Fehler beim Allokieren der Generatoren\n");
		fprintf(stderr,"Fehler beim Allokieren der Generatoren\n");
	}
	#pragma omp parallel for
	for(int core=0;core<maxcores;core+=1){//allokieren und initialisieren aller Generatoren, für jeden Möglichen Kern ein Generator benoetigt
		generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(generatoren[core], seed+core);
	}
	printf("Laenge=%d\tTemperatur=%f\n", laenge, temperatur);
	char gitter[laenge*laenge];//Gitter erstellen und von thermalisieren ausgeben lassen
	initialisierung(gitter, laenge, seed);
	thermalisieren(laenge, temperatur, j, seed, 500, gitter, dummydatei, generatoren[0]);//Erstes Thermalisieren, hier nur zur Ausgabe des Gitters
	//500 Thermalisierungsschritte bei vergleich verschiedener Gitterlaenge, 50000 bei Vergleich verschiedener Temperaturen
	//Zur Beschleunigung auch thermalisierenmehreregeneratoren moeglich
	
	for (int durchlauf=0; durchlauf<durchlaeufe;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden, Unterschiede in der Laufzeit aufgrund von Hitze etc. herauszumitteln
	//Vergleichsmassstab: Messungen bei einem core	
	//	#pragma omp parallel for//Fuer jede Messreihe neu intialisieren, damit immer dasselbe gemessen wird->Zeigt sich, dass dies zu sprunghaften Ergebnissen in der Skalierung führt
	//	for(int core=0;core<maxcores;core+=1){
	//			generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
	//			gsl_rng_set(generatoren[core], seed+core);
	//	}
		speedup=1;//aus definition
		omp_set_num_threads(1);
		einlesen(gitter, laenge, dummydatei);
		messdatei=fopen(dateinamedummymessen,  "w+");
		gettimeofday(&anfangmessen, NULL);//Zeitmessen der Messung
		messenmehreregeneratoren(laenge, temperatur, j, messungen, /*gitter*/dummydatei, messdatei, generatoren);
		gettimeofday(&endemessen, NULL);
		sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
		usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
		zeiteincore=sec+1e-06*usec;//Bestimmen der Zeit, die gebracuht wurde
		ergebnisse[durchlauf]=zeiteincore;
		printf("bei %d cores haben %d Messungen %f Sekunden gebraucht\n", 1, messungen, zeiteincore);
		fprintf(zeitdatei, "%f\t%f\t%f\t%f\t%f\n", 1.0, (double)messungen, zeiteincore, 1.0, (double)laenge);//cores messungen Zeit Speedup
		fclose(messdatei);
	}
	mittelzeit=mittelwertarray(ergebnisse, durchlaeufe);//Mitteln um 
	varianzzeit=varianzarray(ergebnisse, durchlaeufe, mittelzeit);
	zeiteincore=mittelzeit;
	varianzeincore=varianzzeit;
	speedupmittel=1;
	speedupfehler=varianzzeit/mittelzeit;
	zeitmineincore=minarray(ergebnisse, durchlaeufe);
	fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 1.0, (double)laenge, mittelzeit, varianzzeit, speedupmittel, speedupfehler, temperatur, zeitmineincore, 1.0);
		//Messungen bei mehreren cores
	for (int cores=2;cores<=maxcores;cores+=1){
		for (int durchlauf=0; durchlauf<durchlaeufe;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden
			omp_set_num_threads(cores);//Alle moeglichen Kernzahlen werden benutzt
		//	#pragma omp parallel for//Fuer jede Messreihe neu intialisieren, damit immer dasselbe gemessen wird->Zeigt sich, dass dies zu sprunghaften Ergebnissen in der Skalierung führt
		//	for(int core=0;core<maxcores;core+=1){
		//		generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
		//		gsl_rng_set(generatoren[core], seed+core);
		//	}
			einlesen(gitter, laenge, dummydatei);
			messdatei=fopen(dateinamedummymessen, "w+");
			gettimeofday(&anfangmessen, NULL);
			messenmehreregeneratoren(laenge, temperatur, j, messungen, /*gitter*/dummydatei, messdatei, generatoren);
			gettimeofday(&endemessen, NULL);
			sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
			usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
			zeitgesmessen=sec+1e-06*usec;
			ergebnisse[durchlauf]=zeitgesmessen;
			speedup=zeitgesmessen/zeiteincore;
			printf("bei %d cores haben %d Messungen %f Sekunden gebraucht\n", cores, messungen, zeitgesmessen);
			fprintf(zeitdatei, "%f\t%f\t%f\t%f\t%f\n", (double)cores, (double)messungen, zeitgesmessen, speedup, (double)laenge);
			fclose(messdatei);
		}
		//Auswertung: Vergleich mit Ergebnissen bei einem Thread
		mittelzeit=mittelwertarray(ergebnisse, durchlaeufe);
		varianzzeit=varianzarray(ergebnisse, durchlaeufe, mittelzeit);
		speedupmittel=mittelzeit/zeiteincore;
		speedupfehler=sqrt(((varianzzeit/zeiteincore)*(varianzzeit/zeiteincore))+(mittelzeit*varianzeincore/zeiteincore/zeiteincore)*(mittelzeit*varianzeincore/zeiteincore/zeiteincore));
		zeitmin=minarray(ergebnisse, durchlaeufe);
		speedupmin=zeitmin/zeitmineincore;
		fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", (double)cores, (double)laenge, mittelzeit, varianzzeit, speedupmittel, speedupfehler, temperatur, zeitmin, speedupmin);
	}
		//Abschaetzungen overhead
	einlesen(gitter, laenge, dummydatei);
	messdatei=fopen(dateinamedummymessen,  "w+");
	gettimeofday(&anfangmessen, NULL);
	einlesen(gitter, laenge, dummydatei);
	double H=hamiltonian(gitter, laenge, j);
	for (int i=0;i<messungen;i+=1){
	fprintf(messdatei, "%f\t", (double)i);
	fprintf(messdatei, "%f\n", H);
	H=gittersummeohnepar(gitter, laenge)/(double)laenge/(double)laenge;}
	gettimeofday(&endemessen, NULL);
	sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
	usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
	zeiteincore=sec+1e-06*usec;
	fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 0.0, (double)laenge, zeiteincore, 0.0,0.0,0.0, temperatur,0.0,0.0);
	fclose(messdatei);

	
	fclose(dummydatei);
	fclose(zeitdatei);
	fclose(mitteldatei);
	
	for(int core=0;core<maxcores;core+=1){
		gsl_rng_free(generatoren[core]);
	}
	free(ergebnisse);
	free(generatoren);
}
