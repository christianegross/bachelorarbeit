//Christiane, 10.07.20
//Skalierungstest von sweepmpi fuer vorgegebene Prozesszahl mit MPI


#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <mpi.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "mpifunktionen.h"


int main(int argc, char **argv){
	//benoetigte Variablen initialisieren
	MPI_Init(&argc, &argv);//Nur einmal im Pogramm moeglich, ganz am Anfang, damit alle Prozesse alle Variablen haben
	int a,b;
	MPI_Comm_size(MPI_COMM_WORLD, &a);//Fuer Aufteilungen benoetigt
	MPI_Comm_rank(MPI_COMM_WORLD, &b);
	const int anzproz=a, myrank=b;
	int laenge;
	double temperatur;
	char merkmal[50];
	if (argc<3){//Zuweisung ueber Kommandozeile moeglich
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
	double mittelzeit, varianzzeit, mittelzeitges, varianzzeitges;//Fuer Durchschnitte
	double zeitmin, zeitmax, zeitminges, zeitmaxges;//Fuer min und max der einzelenen Messungen
	struct timeval anfangmessen, endemessen;//Zeitmessung mit gettimeofday
	double sec, usec, zeitgesmessen;
	const int rechner=0;//1,2 qbig, 0 vm
	const int durchlaeufe=10;//Anzahl an Zeitmessungen
	double ergebnisselokal[durchlaeufe], ergebnisseglobal[durchlaeufe*anzproz], mittelglobal[5*anzproz];//arrays zum speichern der Zeiten, pro Prozess und fuer alle, mittellokal erst bei erstmaliger Verwendung definiert
	FILE *messdatei, *dummydatei, *mitteldatei;//speichern der Messergenbnisse der einzelnen Messungen, der Zeiten und ein nicht benoetigtes nicht thermalisiertes Gitter
	char dateinamemittel[200], dateinamedummytherm[50], dateinamedummymessen[50];//dateinamen
	sprintf(dateinamedummytherm, "dummytherm%.2d.txt", rechner);
	sprintf(dateinamedummymessen, "dummymessen%.2d.txt", rechner);
	sprintf(dateinamemittel,"Messungen/Zeiten/zmittel-m%.6d-proz%.2d-%.2d%s.txt",messungen, rechner, anzproz, merkmal);
	//In Dateien wird nur von einem Prozess geschrieben, daher nur ein Prozessmit "w" oeffnen, "w" kreiert Datei, deshalb warten, dass Datei sicher existiert, bis andere Prozesse Datei oeffnen	
	if(myrank==0){
		messdatei=fopen(dateinamedummymessen,  "w+");//Speichert einzelne Messergebnisse, ergebnisse nicht benoetigt
		dummydatei=fopen(dateinamedummytherm, "w+");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird	
		mitteldatei=fopen(dateinamemittel, "w+");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0){
		messdatei=fopen(dateinamedummymessen,  "r");
		dummydatei=fopen(dateinamedummytherm, "r");
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
	char gitter[laenge*laenge];//Gitter erstellen und von thermalisieren ausgeben lassen
	if(myrank==0){//Initialisieren nur in einem Prozess, danach auf alle verteilen
		initialisierenhomogen(gitter, laenge);
	}
	MPI_Bcast(gitter, laenge*laenge, MPI_CHAR, 0, MPI_COMM_WORLD);
	thermalisierenmpi(10, laenge, temperatur, j, gitter, messdatei, dummydatei, generatoren);//Erstes Thermalisieren, hier nur zur Ausgabe des Gitters, mit 50000 Schritten, falls bei verschiedenen Temperaturen verglichen werden soll
	fclose(messdatei);
	for (int durchlauf=0; durchlauf<durchlaeufe;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden, Unterschiede in der Laufzeit aufgrund von Hitze etc. herauszumitteln
	//In Dateien wird nur von einem Prozess geschrieben, daher nur ein Prozessmit "w" oeffnen, "w" kreiert Datei, deshalb warten, dass Datei sicher existiert, bis andere Prozesse Datei oeffnen	
		if(myrank==0){
			messdatei=fopen(dateinamedummymessen,  "w+");
			}
		MPI_Barrier(MPI_COMM_WORLD);
		if(myrank!=0){messdatei=fopen(dateinamedummymessen,  "r");}
		gettimeofday(&anfangmessen, NULL);//Zeitmessen der Messung
		messenmpi(messungen, laenge, temperatur, j, gitter, messdatei, generatoren);
		gettimeofday(&endemessen, NULL);
		sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
		usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
		zeitgesmessen=sec+1e-06*usec;//Bestimmen der Zeit, die gebracuht wurde
		ergebnisselokal[durchlauf]=zeitgesmessen;//Einzelen Ergebnisse speichern
		fclose(messdatei);
	}
	MPI_Gather(ergebnisselokal, durchlaeufe, MPI_DOUBLE, ergebnisseglobal, durchlaeufe, MPI_DOUBLE,0, MPI_COMM_WORLD);//Alle einzelmessungen sammeln
	mittelzeit=mittelwertarray(ergebnisselokal, durchlaeufe);//Durchschnitt bei einem Prozess
	varianzzeit=varianzarray(ergebnisselokal, durchlaeufe, mittelzeit);//Standardabweichung bei einem Prozess
	zeitmin=minarray(ergebnisselokal, durchlaeufe);//Minimum bei einem Prozess
	zeitmax=maxarray(ergebnisselokal, durchlaeufe);//Maximum bei einem Prozess
	mittelzeitges=mittelwertarray(ergebnisseglobal, durchlaeufe*anzproz);//Durchschnitt ueber alle Prozesse
	varianzzeitges=varianzarray(ergebnisseglobal, durchlaeufe*anzproz, mittelzeitges);//Standardabweichung ueber alle Prozesse
	zeitminges=minarray(ergebnisseglobal, durchlaeufe*anzproz);//Minimum aller Prozesse
	zeitmaxges=maxarray(ergebnisseglobal, durchlaeufe*anzproz);//Maximum aller Prozesse
	double mittellokal[5]={(double)myrank, mittelzeit, varianzzeit, zeitmin, zeitmax};
	MPI_Gather(mittellokal, 5, MPI_DOUBLE, mittelglobal, 5, MPI_DOUBLE,0 , MPI_COMM_WORLD);//Ergebnisse von allen Prozessen sammeln
	if(myrank==0){
		for(int i=0; i<anzproz; i+=1){//Ergebnisse der einzelnen Prozesse ausgeben
			fprintf(mitteldatei,  "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",
				mittelglobal[5*i], (double)anzproz, (double)temperatur, (double)laenge, mittelglobal[5*i+1], mittelglobal[5*i+2], mittelglobal[5*i+3], mittelglobal[5*i+4]); 
		}
		fprintf(mitteldatei,  "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",//Ergebniss aller Prozesse ausgeben
				(double)anzproz, (double)anzproz, (double)temperatur, (double)laenge, mittelzeitges, varianzzeitges, zeitminges, zeitmaxges); 
	}	
	ausgabe(gitter, laenge, dummydatei);
	fclose(dummydatei);
	fclose(mitteldatei);
	for(int core=0;core<anzproz;core+=1){
		gsl_rng_free(generatoren[core]);
	}
	free(generatoren);
	MPI_Finalize();
}
