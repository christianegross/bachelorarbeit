//Christiane, 25.05.20
//Skalierungstest von sweep fuer verschiedene cores


#include <stdio.h>
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
	int laenge;//laenge der verwendeten Gitter
	int lenarray[3]={10, 100, 500};//, 10, 50, 50};
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int messungen=1000;//pro temperatur
	double mittelzeit, varianzzeit, speedupmittel, speedupfehler, speedup;
	int node=2;//1,2 qbig, 0 vm
	char merkmal[50]="mg2snichtimmerzufall-mehreret";
	int durchlaeufe=5;
	double temperatur=4.5;//Skalierung bei nur einer Temperatur messen niedrig 0.5, mittel2, mittel2 2.5, hoch 3.5
	double temperaturen[4]={0.5, 2, 2.5, 4.5};
	//t1=0.5, t2=2, t3=2.5, t4=4.5
	double *ergebnisse;
	if((ergebnisse=(double*)malloc(sizeof(double)*durchlaeufe))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren der Temperaturen!\n");
		return (-1);
	}
	//~ double mittelham, varianzham;
	FILE *messdatei, *dummydatei/*, *dummydateiplot*/, *zeitdatei, *mitteldatei/*, *vergleichsdatei*/;//speichern der Messergenbnisse der einzelnen Messungen, der Zeiten und ein nicht benoetigtes nicht thermalisiertes Gitter
	char dateinamezeit[200], dateinamemittel[200], dateinamedummytherm[50]/*, dateinamedummythermplot[50]*/, dateinamedummymessen[50];//dateinamemessen[150],
	struct timeval anfangmessen, endemessen;//Zeitmessung mit gettimeofday
	double sec, usec, zeitgesmessen, zeiteincore, varianzeincore;
	sprintf(dateinamedummytherm, "dummytherm%.2d.txt", node);
	//sprintf(dateinamedummythermplot, "dummythermplot%.2d.txt", node);
	sprintf(dateinamedummymessen, "dummymessen%.2d.txt", node);

	dummydatei=fopen(dateinamedummytherm, "w+");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird	
	//dummydateiplot=fopen(dateinamedummythermplot, "w+");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird	
	//sprintf(dateinamezeit,"Messungen/Zeiten/zeitenmessen-laenge%.4d-m%.6d-mehrere.txt",laenge,messungen);
	sprintf(dateinamezeit,"Messungen/Zeiten/zmessen-m%.6d-node%.2d%s.txt",messungen, node, merkmal);
	sprintf(dateinamemittel,"Messungen/Zeiten/zmittel-m%.6d-node%.2d%s.txt",messungen, node, merkmal);
	zeitdatei=fopen(dateinamezeit, "w+");
	mitteldatei=fopen(dateinamemittel, "w+");
	//gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	//gsl_rng_set(generator, seed);
	//Versuch: jeder thread einen einzelnen Generator
	gsl_rng **generatoren;
	if ((generatoren=(gsl_rng**)malloc(maxcores*sizeof(gsl_rng**)))==NULL){
		printf("Fehelr beim Allokieren der Generatoren\n");
	}
	#pragma omp parallel for
	for(int core=0;core<maxcores;core+=1){
		generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(generatoren[core], seed+core);
	}
	for(int t=0;t<4;t+=1){
		temperatur=temperaturen[t];
	for (int laengen=0; laengen<3; laengen+=1){
		laenge=lenarray[laengen];
		printf("Laenge=%d\n", laenge);
		char gitter[laenge*laenge];//Gitter erstellen und von thermalisieren ausgeben lassen
		initialisierung(gitter, laenge, seed);
		thermalisieren(laenge, temperatur, j, seed, 500, gitter, dummydatei, generatoren[0]);//Erstes Thermalisieren, hier nur zur Ausgabe des Gitters
		for (int durchlauf=0; durchlauf<durchlaeufe;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden
		//Vergleichsmassstab: Messungen bei einem core	
			//gsl_rng_set(generator, seed/*(seed+durchlauf)*/);
			#pragma omp parallel for
			for(int core=0;core<maxcores;core+=1){
					generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
					gsl_rng_set(generatoren[core], seed+core);
			}
			speedup=1;//aus definition
			omp_set_num_threads(1);
			einlesen(gitter, laenge, dummydatei);
			//sprintf(dateinamemessen,"Messungen/Messwerte/messung-laenge%.4d-m%.6d-cores%.2d-%.2d.txt",laenge,messungen,1, durchlauf);
			//messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
			messdatei=fopen(dateinamedummymessen,  "w+");
			//~ vergleichsdatei=fopen("dummyvergleich.txt", "w+");
			gettimeofday(&anfangmessen, NULL);
			//~ messenvergleichen(laenge, temperatur, j, messungen, dummydatei, messdatei, vergleichsdatei, generator);
			messenmehreregeneratoren(laenge, temperatur, j, messungen, gitter/*dummydatei*/, messdatei, generatoren);
			gettimeofday(&endemessen, NULL);
			//~ mittelham=mittelwertberechnungnaiv(vergleichsdatei, messungen, 3, 4);
			//~ varianzham=varianzberechnungnaiv(vergleichsdatei, messungen, mittelham, 3, 4);
			//~ printf("Laenge %d, DeltaH=%f+/-%f\n", laenge, mittelham, varianzham);
			sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
			usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
			zeiteincore=sec+1e-06*usec;
			ergebnisse[durchlauf]=zeiteincore;
			printf("bei %d cores haben %d Messungen %f Sekunden gebraucht\n", 1, messungen, zeiteincore);
			fprintf(zeitdatei, "%f\t%f\t%f\t%f\t%f\n", 1.0, (double)messungen, zeiteincore, 1.0, (double)laenge);//cores messungen Zeit Speedup
			fclose(messdatei);
			//~ fclose(vergleichsdatei);
		}
		mittelzeit=mittelwertarray(ergebnisse, durchlaeufe);
		varianzzeit=varianzarray(ergebnisse, durchlaeufe, mittelzeit);
		zeiteincore=mittelzeit;
		varianzeincore=varianzzeit;
		speedupmittel=1;
		speedupfehler=varianzzeit/mittelzeit;
		fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 1.0, (double)laenge, mittelzeit, varianzzeit, speedupmittel, speedupfehler, temperatur);
			//Messungen bei mehreren cores
		for (int cores=2;cores<=maxcores;cores+=1){
			for (int durchlauf=0; durchlauf<durchlaeufe;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden
				omp_set_num_threads(cores);
				einlesen(gitter, laenge, dummydatei);
				//sprintf(dateinamemessen,"Messungen/Messwerte/messung-laenge%.4d-m%.6d-cores%.2d-%.2d.txt",laenge,messungen,cores,durchlauf);
				//messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
				messdatei=fopen(dateinamedummymessen, "w+");
				//~ vergleichsdatei=fopen("dummyvergleich.txt", "w+");
				gettimeofday(&anfangmessen, NULL);
				//~ messenvergleichen(laenge, temperatur, j, messungen, dummydatei, messdatei, vergleichsdatei, generator);
				messenmehreregeneratoren(laenge, temperatur, j, messungen, gitter/*dummydatei*/, messdatei, generatoren);
				gettimeofday(&endemessen, NULL);
				//~ mittelham=mittelwertberechnungnaiv(vergleichsdatei, messungen, 3, 4);
				//~ varianzham=varianzberechnungnaiv(vergleichsdatei, messungen, mittelham, 3, 4);
				//~ printf("Laenge %d, DeltaH=%f+/-%f\n", laenge, mittelham, varianzham);
				sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
				usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
				zeitgesmessen=sec+1e-06*usec;
				ergebnisse[durchlauf]=zeitgesmessen;
				speedup=zeitgesmessen/zeiteincore;
				printf("bei %d cores haben %d Messungen %f Sekunden gebraucht\n", cores, messungen, zeitgesmessen);
				fprintf(zeitdatei, "%f\t%f\t%f\t%f\t%f\n", (double)cores, (double)messungen, zeitgesmessen, speedup, (double)laenge);
				fclose(messdatei);
		//~ fclose(vergleichsdatei);
			}
			mittelzeit=mittelwertarray(ergebnisse, durchlaeufe);
			varianzzeit=varianzarray(ergebnisse, durchlaeufe, mittelzeit);
			speedupmittel=mittelzeit/zeiteincore;
			speedupfehler=sqrt(((varianzzeit/zeiteincore)*(varianzzeit/zeiteincore))+(mittelzeit*varianzeincore/zeiteincore/zeiteincore)*(mittelzeit*varianzeincore/zeiteincore/zeiteincore));
			fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", (double)cores, (double)laenge, mittelzeit, varianzzeit, speedupmittel, speedupfehler, temperatur);
		}
			//Messungen overhead
		einlesen(gitter, laenge, dummydatei);
		messdatei=fopen(dateinamedummymessen,  "w+");
		gettimeofday(&anfangmessen, NULL);
		//einlesen(gitter, laenge, dummydatei);
		double H=hamiltonian(gitter, laenge, j);
		for (int i=0;i<messungen;i+=1){
		fprintf(messdatei, "%f\t", (double)i);
		fprintf(messdatei, "%f\n", H);
		H=gittersummeohnepar(gitter, laenge)/(double)laenge/(double)laenge;}
		gettimeofday(&endemessen, NULL);
		sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
		usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
		zeiteincore=sec+1e-06*usec;
		fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 0.0, (double)laenge, zeiteincore, 0.0,0.0,0.0, temperatur);
		fclose(messdatei);
	}
	}

	
	fclose(dummydatei);
	fclose(zeitdatei);
	fclose(mitteldatei);
	
	//gsl_rng_free(generator);//free, close: zum Verhindern von Speicherproblemen
	for(int core=0;core<maxcores;core+=1){
		gsl_rng_free(generatoren[core]);
	}
	free(ergebnisse);
	free(generatoren);
}
