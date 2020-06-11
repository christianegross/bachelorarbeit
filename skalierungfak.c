//Ausprobieren der Skalierung auf QBiG
//Christiane, 11.06.20

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>//Parallelisierung
#include "math.h"
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen

int fak(int argument){
	int ergebnis=1;
	if(argument<0){
		printf("Ungueltiges Argument!\n");
		ergebnis=-1;
	}
	else{
		for (int i=1; i<=argument; i+=1){
			ergebnis*=i;
		}
	}
	return ergebnis;
}

double mittelwertarray(double *array, int messungen){
	//Bestimmt Mittelwert eines arrays, das mit doubles gefüllt ist
	double summe=0;//Speichert Summe über 
	for (int messung=0; messung<messungen; messung+=1){//Mittelwert über Messung bilden
		summe+=array[messung];
	}
	return summe/(double)messungen;
}

double varianzarray(double *array, int messungen, double mittelwert){
	//Bestimmt die Varianz über alle Elemente in einem array, das mit doubles gefuellt ist, um den Mittelwert
	double summe=0;//Speichert Summe über Messungen
	double einwert=0;//Speichert einen ausgelesenen Wert
	for (int messung=0; messung<messungen; messung+=1){//Mittelwert über Messung bilden
		einwert=array[messung];
		summe+=(einwert-mittelwert)*(einwert-mittelwert);
	}
	return sqrt(summe/((double)messungen-1));
}


int main(int argc, char **argv){
	int maxcores=omp_get_max_threads();//aus Computerarchitektut/batchskript
	int messungen=100000;
	int argument=20000;
	int node=1;//1,2 qbig, 0 vm
	int durchlaeufe=10;
	double mittelzeit, varianzzeit;
	double *ergebnisse;
	if((ergebnisse=(double*)malloc(sizeof(double)*durchlaeufe))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren der Temperaturen!\n");
		return (-1);
	}
	FILE *zeitdatei;
	char dateinamezeit[200];
	struct timeval anfangmessen, endemessen;//Zeitmessung mit gettimeofday
	double sec, usec, zeitgesmessen;
	sprintf(dateinamezeit,"Messungen/Zeiten/skalierungfak-m-%.8d-node%.2d",messungen, node);
	zeitdatei=fopen(dateinamezeit, "w+");
	int ergebnis=fak(argument);
	for (int cores=1;cores<=maxcores;cores+=1){
			for (int durchlauf=0; durchlauf<durchlaeufe;durchlauf+=1){//mehrere Durchläufe, um Unstimmigkeiten mit gettimeofday herauszufinden
				printf("%d\n", durchlauf);
				omp_set_num_threads(cores);
				gettimeofday(&anfangmessen, NULL);
				#pragma omp parallel for
				for (int i=0; i<messungen;i+=1){
					if(ergebnis!=fak(argument)){
						printf("Fehler!\n");
					}
				}
				gettimeofday(&endemessen, NULL);
				sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
				usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
				zeitgesmessen=sec+1e-06*usec;
				ergebnisse[durchlauf]=zeitgesmessen;
			}
			mittelzeit=mittelwertarray(ergebnisse, durchlaeufe);
			varianzzeit=varianzarray(ergebnisse, durchlaeufe, mittelzeit);
			fprintf(zeitdatei, "%f\t%f\t%f\t%f\n", (double)argument, (double)cores, mittelzeit, varianzzeit);
		}
	fclose(zeitdatei);
	free(ergebnisse);
}
