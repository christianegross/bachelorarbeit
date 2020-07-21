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
	int maxproz, myrank;
	MPI_Comm_size(MPI_COMM_WORLD, &maxproz);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
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
	double mittelzeit, varianzzeit, mittelzeitges, varianzzeitges;
	double zeiteinproz, varianzeinproz, speedup, speedupvar;
	double zeitmin, zeitmax, zeitminges, zeitmaxges;
	struct timeval anfangmessen, endemessen;//Zeitmessung mit gettimeofday
	double sec, usec, zeitgesmessen;
	const int rechner=0;//1,2 qbig, 0 vm
	const int durchlaeufe=10;
	double ergebnisselokal[durchlaeufe];//, ergebnisseglobal[durchlaeufe*anzproz], mittelglobal[5*anzproz];//mittellokal erst bei erstmaliger Verwendung definiert
	double *ergebnisseglobal, *mittelglobal, mittellokal[5]={1,1,1,1,1};
	if((ergebnisseglobal=(double*)malloc(sizeof(double)*durchlaeufe*maxproz))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren von ergebnisseglobal!\n");
		fprintf(stderr, "Fehler beim Allokieren ergebnisseglobal!\n");
		return (-1);
	}
	if((mittelglobal=(double*)malloc(sizeof(double)*5*maxproz))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren von mittelglobal!\n");
		fprintf(stderr, "Fehler beim Allokieren mittelglobal!\n");
		return (-1);
	}
	FILE *messdatei, *dummydatei, *mitteldatei, *mitteldateiglobal;;//speichern der Messergenbnisse der einzelnen Messungen, der Zeiten und ein nicht benoetigtes nicht thermalisiertes Gitter
	char dateinamemittelglobal[200], dateinamemittel[200], dateinamedummytherm[50], dateinamedummythermplot[50], dateinamedummymessen[50];//dateinamemessen[150],
	sprintf(dateinamedummytherm, "dummytherm%.2d.txt", rechner);
	sprintf(dateinamedummymessen, "dummymessen%.2d.txt", rechner);
	//sprintf(dateinamezeit,"Messungen/Zeiten/zmessen-m%.6d-proz%.2d-%.2d%s.txt",messungen, rechner, anzproz, merkmal);
	sprintf(dateinamemittel,"Messungen/Zeiten/zmittel-m%.6d-proz%.2d-%.2d%s.txt",messungen, rechner, maxproz, merkmal);
	sprintf(dateinamemittelglobal,"Messungen/Zeiten/zmittelglobal-m%.6d-proz%.2d-%.2d%s.txt",messungen, rechner, maxproz, merkmal);
	if(myrank==0){
		messdatei=fopen(dateinamedummymessen,  "w+");
		dummydatei=fopen(dateinamedummytherm, "w+");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird	
		//zeitdatei=fopen(dateinamezeit, "w+");
		mitteldatei=fopen(dateinamemittel, "w+");
		mitteldateiglobal=fopen(dateinamemittelglobal, "w+");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0){
		messdatei=fopen(dateinamedummymessen,  "r");
		dummydatei=fopen(dateinamedummytherm, "r");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird	
		//zeitdatei=fopen(dateinamezeit, "r");
		mitteldatei=fopen(dateinamemittel, "r");
		mitteldateiglobal=fopen(dateinamemittelglobal, "r");
	}
	gsl_rng **generatoren;
	if ((generatoren=(gsl_rng**)malloc(maxproz*sizeof(gsl_rng**)))==NULL){
		printf("Fehler beim Allokieren der Generatoren\n");
		fprintf(stderr,"Fehler beim Allokieren der Generatoren\n");
	}
	for(int core=0;core<maxproz;core+=1){//allokieren und initialisieren aller Generatoren
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
	//printf("3 %d\n", myrank);
	thermalisierenmpi(10, laenge, maxproz, myrank, temperatur, j, gitter, messdatei, dummydatei, generatoren);//Erstes Thermalisieren, hier nur zur Ausgabe des Gitters
	//printf("4 Laenge=%d\tTemperatur=%f\n", laenge, temperatur);
	fclose(messdatei);
	//Messungen mit einem Prozess
	if (myrank==0){
		for (int durchlauf=0; durchlauf<durchlaeufe;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden
			//einlesen(gitter, laenge, dummydatei);
			if(myrank==0){
				messdatei=fopen(dateinamedummymessen,  "w+");
				//printf("%d\n", durchlauf);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if(myrank!=0){messdatei=fopen(dateinamedummymessen,  "r");}
			
		//printf("5 %d\n", myrank);
			gettimeofday(&anfangmessen, NULL);
			messenmpi(messungen, laenge, 1, myrank, temperatur, j, gitter, messdatei, generatoren);
		//printf("6 %d\n", myrank);
			gettimeofday(&endemessen, NULL);
			sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
			usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
			zeitgesmessen=sec+1e-06*usec;
			ergebnisselokal[durchlauf]=zeitgesmessen;
			//printf("bei core %d von %d haben %d Messungen %f Sekunden gebraucht\n", myrank, anzproz, messungen, zeiteincore);
			printf("%f\t%f\t%f\t%f\t%f\n", (double)myrank, 1.0, (double)messungen, (double)laenge, zeitgesmessen);//cores messungen Zeit Speedup
			fclose(messdatei);
		}
		//printf("vor gather2 in %d\n", myrank);
		mittelzeit=mittelwertarray(ergebnisselokal, durchlaeufe);
		varianzzeit=varianzarray(ergebnisselokal, durchlaeufe, mittelzeit);
		zeitmin=minarray(ergebnisselokal, durchlaeufe);
		zeitmax=maxarray(ergebnisselokal, durchlaeufe);
		mittellokal[0]=(double)myrank;
		mittellokal[1]=mittelzeit;
		mittellokal[2]=varianzzeit;
		mittellokal[3]=zeitmin;
		mittellokal[4]=zeitmax;
	}
	MPI_Gather(ergebnisselokal, durchlaeufe, MPI_DOUBLE, ergebnisseglobal, durchlaeufe, MPI_DOUBLE,0, MPI_COMM_WORLD);
	MPI_Gather(mittellokal, 5, MPI_DOUBLE, mittelglobal, 5, MPI_DOUBLE,0 , MPI_COMM_WORLD);
	//printf("nach gather in %d\n", myrank);
	if(myrank==0){
		mittelzeitges=mittelwertarray(ergebnisseglobal, durchlaeufe*1);
		varianzzeitges=varianzarray(ergebnisseglobal, durchlaeufe*1, mittelzeitges);
		zeitminges=minarray(ergebnisseglobal, durchlaeufe*1);
		zeitmaxges=maxarray(ergebnisseglobal, durchlaeufe*1);
		zeiteinproz=mittelzeitges;
		varianzeinproz=varianzzeitges;
		for(int i=0; i<1; i+=1){
			fprintf(mitteldatei,  "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
				mittelglobal[5*i], 1.0, (double)temperatur, (double)laenge, mittelglobal[5*i+1], mittelglobal[5*i+2], mittelglobal[5*i+3], mittelglobal[5*i+4]); 
		}
		fprintf(mitteldateiglobal,  "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
				(double)maxproz, (double)temperatur, (double)laenge, mittelzeitges, varianzzeitges, zeitminges, zeitmaxges, 1.0, varianzzeitges/mittelzeitges); 
	}	
	for(int anzproz=2;anzprozymaxproz;anzproz+=1){
		for (int durchlauf=0; durchlauf<durchlaeufe;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden
			//einlesen(gitter, laenge, dummydatei);
			if(myrank==0){
				messdatei=fopen(dateinamedummymessen,  "w+");
				//printf("%d\n", durchlauf);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if(myrank!=0){messdatei=fopen(dateinamedummymessen,  "r");}
			
		//printf("5 %d\n", myrank);
			gettimeofday(&anfangmessen, NULL);
			messenmpi(messungen, laenge, 1, myrank, temperatur, j, gitter, messdatei, generatoren);
		//printf("6 %d\n", myrank);
			gettimeofday(&endemessen, NULL);
			sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
			usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
			zeitgesmessen=sec+1e-06*usec;
			ergebnisselokal[durchlauf]=zeitgesmessen;
			//printf("bei core %d von %d haben %d Messungen %f Sekunden gebraucht\n", myrank, anzproz, messungen, zeiteincore);
			printf("%f\t%f\t%f\t%f\t%f\n", (double)myrank, 1.0, (double)messungen, (double)laenge, zeitgesmessen);//cores messungen Zeit Speedup
			fclose(messdatei);
		}
		//printf("vor gather2 in %d\n", myrank);
		mittelzeit=mittelwertarray(ergebnisselokal, durchlaeufe);
		varianzzeit=varianzarray(ergebnisselokal, durchlaeufe, mittelzeit);
		zeitmin=minarray(ergebnisselokal, durchlaeufe);
		zeitmax=maxarray(ergebnisselokal, durchlaeufe);
		mittellokal[0]=(double)myrank;
		mittellokal[1]=mittelzeit;
		mittellokal[2]=varianzzeit;
		mittellokal[3]=zeitmin;
		mittellokal[4]=zeitmax;
	}
	MPI_Gather(ergebnisselokal, durchlaeufe, MPI_DOUBLE, ergebnisseglobal, durchlaeufe, MPI_DOUBLE,0, MPI_COMM_WORLD);
	MPI_Gather(mittellokal, 5, MPI_DOUBLE, mittelglobal, 5, MPI_DOUBLE,0 , MPI_COMM_WORLD);
	//printf("nach gather in %d\n", myrank);
	if(myrank==0){
		mittelzeitges=mittelwertarray(ergebnisseglobal, durchlaeufe*1);
		varianzzeitges=varianzarray(ergebnisseglobal, durchlaeufe*1, mittelzeitges);
		zeitminges=minarray(ergebnisseglobal, durchlaeufe*1);
		zeitmaxges=maxarray(ergebnisseglobal, durchlaeufe*1);
		zeiteinproz=mittelzeitges;
		varianzeinproz=varianzzeitges;
		for(int i=0; i<1; i+=1){
			fprintf(mitteldatei,  "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
				mittelglobal[5*i], 1.0, (double)temperatur, (double)laenge, mittelglobal[5*i+1], mittelglobal[5*i+2], mittelglobal[5*i+3], mittelglobal[5*i+4]); 
		}
		fprintf(mitteldateiglobal,  "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
				(double)maxproz, (double)temperatur, (double)laenge, mittelzeitges, varianzzeitges, zeitminges, zeitmaxges, 1.0, varianzzeitges/mittelzeitges); 
	}	
	//printf("Ausgabe in %d\n", myrank);
	ausgabe(gitter, laenge, dummydatei);
	fclose(dummydatei);
	//fclose(zeitdatei);
	fclose(mitteldatei);
	fclose(mitteldateiglobal);
	for(int core=0;core<maxproz;core+=1){
		gsl_rng_free(generatoren[core]);
	}
	free(generatoren);
	free(mittelglobal);
	free(ergebnisseglobal);
	MPI_Finalize();
}
