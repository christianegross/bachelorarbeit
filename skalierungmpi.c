//Christiane, 10.07.hproz20
//Skalierungstest von sweep fuer verschiedene cores mit MPI


#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <mpi.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "mpifunktionen.h"


int main(int argc, char **argv){
	//benoetigte Variablen initialisieren
	//int maxcores=omp_get_max_threads();//aus Computerarchitektut/batchskript
	MPI_Init(&argc, &argv);//Anzahl Prozessoren, meine Nummer
	int a,b;
	MPI_Comm_size(MPI_COMM_WORLD, &a);
	MPI_Comm_rank(MPI_COMM_WORLD, &b);
	const int anzproz=a, myrank=b;
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
	else{	
	laenge=atoi(argv[1]);//laenge der verwendeten Gitter
	temperatur=atof(argv[2]);//Temperatur
	sprintf(merkmal,"%s-l%s",argv[3], argv[1]);
	}
	double j=1.0;
	const int seed=5;//fuer den zufallsgenerator
	const int messungen=1000;//pro temperatur
	double mittelzeit, varianzzeit;
	double zeitmin, zeitmax;
	struct timeval anfangmessen, endemessen;//Zeitmessung mit gettimeofday
	double sec, usec, zeitgesmessen;
	const int rechner=0;//1,2 qbig, 0 vm
	const int durchlaeufe=10;
	//~ double *ergebnisselokal, *ergebnisseglobal, *mittelglobal;//mittellokal erst bei erstmaliger Verwendung definiert
	double ergebnisselokal[durchlaeufe], ergebnisseglobal[durchlaeufe*anzproz], mittelglobal[5*anzproz];//mittellokal erst bei erstmaliger Verwendung definiert
	//~ if((ergebnisselokal=(double*)malloc(sizeof(double)*durchlaeufe))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		//~ printf("Fehler beim Allokieren der Temperaturen!\n");
		//~ fprintf(stderr, "Fehler beim Allokieren der Temperaturen!\n");
		//~ return (-1);
	//~ }
	//~ if((ergebnisseglobal=(double*)malloc(sizeof(double)*durchlaeufe*anzproz))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		//~ printf("Fehler beim Allokieren der Temperaturen!\n");
		//~ fprintf(stderr, "Fehler beim Allokieren der Temperaturen!\n");
		//~ return (-1);
	//~ }
	//~ if((mittelglobal=(double*)malloc(sizeof(double)*5*anzproz))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		//~ printf("Fehler beim Allokieren der Temperaturen!\n");
		//~ fprintf(stderr, "Fehler beim Allokieren der Temperaturen!\n");
		//~ return (-1);
	//~ }
	FILE *messdatei, *dummydatei, *zeitdatei, *mitteldatei;//speichern der Messergenbnisse der einzelnen Messungen, der Zeiten und ein nicht benoetigtes nicht thermalisiertes Gitter
	char dateinamezeit[200], dateinamemittel[200], dateinamedummytherm[50]/*, dateinamedummythermplot[50]*/, dateinamedummymessen[50];//dateinamemessen[150],
	sprintf(dateinamedummytherm, "dummytherm%.2d.txt", rechner);
	sprintf(dateinamedummymessen, "dummymessen%.2d.txt", rechner);
	sprintf(dateinamezeit,"Messungen/Zeiten/zmessen-m%.6d-proz%.2d-%.2d%s.txt",messungen, rechner, anzproz, merkmal);
	sprintf(dateinamemittel,"Messungen/Zeiten/zmittel-m%.6d-proz%.2d-%.2d%s.txt",messungen, rechner, anzproz, merkmal);
	if(myrank==0){
		messdatei=fopen(dateinamedummymessen,  "w+");
		dummydatei=fopen(dateinamedummytherm, "w+");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird	
		zeitdatei=fopen(dateinamezeit, "w+");
		mitteldatei=fopen(dateinamemittel, "w+");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0){
		messdatei=fopen(dateinamedummymessen,  "r");
		dummydatei=fopen(dateinamedummytherm, "r");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird	
		zeitdatei=fopen(dateinamezeit, "r");
		mitteldatei=fopen(dateinamemittel, "r");
	}
	gsl_rng **generatoren;
	if ((generatoren=(gsl_rng**)malloc(anzproz*sizeof(gsl_rng**)))==NULL){
		printf("Fehler beim Allokieren der Generatoren\n");
		fprintf(stderr,"Fehler beim Allokieren der Generatoren\n");
	}
	for(int core=0;core<anzproz;core+=1){//allokieren und initialisieren aller Generatoren
		generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(generatoren[core], seed+core);
	}
	//printf("1 Laenge=%d\tTemperatur=%f\n", laenge, temperatur);
	int gitter[laenge*laenge];//Gitter erstellen und von thermalisieren ausgeben lassen
	//printf("2 Laenge=%d\tTemperatur=%f\n", laenge, temperatur);	
	if(myrank==0){
		initialisierenhomogen(gitter, laenge);
	}
	MPI_Bcast(gitter, laenge*laenge, MPI_INT, 0, MPI_COMM_WORLD);
	thermalisierenmpi(messungen, laenge, temperatur, j, gitter, messdatei, dummydatei, generatoren);//Erstes Thermalisieren, hier nur zur Ausgabe des Gitters
	//~ printf("4 Laenge=%d\tTemperatur=%f\n", laenge, temperatur);
	for (int durchlauf=0; durchlauf<durchlaeufe;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden
		//einlesen(gitter, laenge, dummydatei);
		if(myrank==0){messdatei=fopen(dateinamedummymessen,  "w+");}
		MPI_Barrier(MPI_COMM_WORLD);
		if(myrank!=0){messdatei=fopen(dateinamedummymessen,  "r");}
		gettimeofday(&anfangmessen, NULL);
		messenmpi(messungen, laenge, temperatur, j, gitter, messdatei, generatoren);
		gettimeofday(&endemessen, NULL);
		sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
		usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
		zeitgesmessen=sec+1e-06*usec;
		ergebnisselokal[durchlauf]=zeitgesmessen;
		//printf("bei core %d von %d haben %d Messungen %f Sekunden gebraucht\n", myrank, anzproz, messungen, zeiteincore);
		printf("%f\t%f\t%f\t%f\t%f\n", (double)myrank, (double)anzproz, (double)messungen, (double)laenge, zeitgesmessen);//cores messungen Zeit Speedup
		fclose(messdatei);
	}
	printf("vor gather2 in %d\n", myrank);
	MPI_Allgather(ergebnisselokal, durchlaeufe, MPI_DOUBLE, ergebnisseglobal, durchlaeufe*anzproz, MPI_DOUBLE, MPI_COMM_WORLD);
	mittelzeit=mittelwertarray(ergebnisselokal, durchlaeufe);
	varianzzeit=varianzarray(ergebnisselokal, durchlaeufe, mittelzeit);
	zeitmin=minarray(ergebnisselokal, durchlaeufe);
	zeitmax=maxarray(ergebnisselokal, durchlaeufe);
	double mittellokal[5]={(double)myrank, mittelzeit, varianzzeit, zeitmin, zeitmax};
	printf("vor gather in %d\n", myrank);
	MPI_Allgather(mittellokal, 5, MPI_DOUBLE, mittelglobal, 5*anzproz, MPI_DOUBLE, MPI_COMM_WORLD);
	
	printf("nach gather in %d\n", myrank);
	if(myrank==0){
		for(int i=0; i<anzproz; i+=1){
			fprintf(mitteldatei,  "%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
				mittelglobal[5*i], (double)anzproz, (double)temperatur, mittelglobal[5*i+1], mittelglobal[5*i+2], mittelglobal[5*i+3], mittelglobal[5*i+4]); 
		}
	}
	printf("Ausgabe in %d\n", myrank);
	//fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\n", (double)laenge, temperatur, mittelzeit, varianzzeit, zeitmineincore);
	
		//Messungen overhead
	//einlesen(gitter, laenge, dummydatei);
	//~ if(myrank==0){messdatei=fopen(dateinamedummymessen,  "w+");}
	//~ MPI_Barrier(MPI_COMM_WORLD);
	//~ if(myrank!=0){messdatei=fopen(dateinamedummymessen,  "r");}
	//~ gettimeofday(&anfangmessen, NULL);
	//~ einlesen(gitter, laenge, dummydatei);
	//~ double H=hamiltonian(gitter, laenge, j);
	//~ for (int i=0;i<messungen;i+=1){
		//~ fprintf(messdatei, "%f\t", (double)i);
		//~ fprintf(messdatei, "%f\n", H);
		//~ H=hamiltonian(gitter, laenge, j);
	//~ }
	//~ gettimeofday(&endemessen, NULL);
	//~ sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
	//~ usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
	//~ zeiteincore=sec+1e-06*usec;
	//~ fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", 0.0, (double)laenge, zeiteincore, 0.0,0.0,0.0, temperatur,0.0,0.0);
	//~ fclose(messdatei);

	
	printf("vor close in %d\n", myrank);
	fclose(dummydatei);
	fclose(zeitdatei);
	fclose(mitteldatei);
	
	for(int core=0;core<anzproz;core+=1){
		gsl_rng_free(generatoren[core]);
	}
	
	printf("vor free malloc in %d\n", myrank);
	free(generatoren);
	MPI_Finalize();
}
