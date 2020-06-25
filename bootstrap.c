//Christiane, start 15.04.20
//Erste Implementierung Ising-Modell für Bachelorarbeit
//Versuch Metropolis-Algorithmus


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
	int anzahlcores=12;
	omp_set_num_threads(anzahlcores);//Setzt die nummer an Kernen, die in den parallelen Regionen verwendet werden.
	int laenge=120;//laenge der verwendeten Gitter
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int N01=10000;//sweeps beim ersten Thermalisieren
	int N0=5000;//benoetigte sweeps zum Thermalisieren
	int messungen=10000;//pro temperatur, zweierpotenz um blocken einfacher zu machen
	int r;//Anzahl an samples für den Bootstrap
	FILE *messdatei, *bootstrapalledatei, *ableitungdatei, *zeitdatei;//benoetigte Dateien zur Ausgabe
	int temperaturzahl=500;//Temperaturen, beid enen gemessen wird
	int schritt=1;//Wie viele Punkte werden gemessen?
	int node=2;//nodes auf vm, qbig
	char dateinamemessen[150], dateinamebootstrapalle[150], dateinameableitung[150], dateinamezeit[150];//Um Dateien mit Variablen benennen zu koennen
	double *temperaturarray;
	if((temperaturarray=(double*)malloc(sizeof(double)*temperaturzahl))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren der Temperaturen!\n");
		return (-1);
	}
	for (int i=0; i<temperaturzahl;i++){//Temperaturarray intalisieren
		temperaturarray[i]=0.01*i+0.01;
		if (i>250){temperaturarray[i]+=(i-250)*0.01;}
	}
	int l;//Laenge der Blocks
	double *blockarray;//Zum Speichern der geblockten Messwerte
	double blocklenarray[12]={32, 64,128, 256, 384, 512, 640, 758, 876, 1024, 1280, 1536};//Blocklaengen, bei denen gemessen wird
	//gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	gsl_rng **generatoren;//mehrere generaotren fuer parallelisierung
	if ((generatoren=(gsl_rng**)malloc(anzahlcores*sizeof(gsl_rng**)))==NULL){
		printf("Fehler beim Allokieren der Generatoren\n");
	}
	for(int core=0;core<anzahlcores;core+=1){
		generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(generatoren[core], seed+core);
	}
	sprintf(dateinamebootstrapalle,"Messungen/Bootstrapges/bootstrapalle-l%.4d-m-%.6d-node%.2d.txt",laenge, messungen, node);//speichert Mitteelwerte aus Bootstrap
	sprintf(dateinameableitung,"Messungen/ableitung-laenge-%.4d-m-%.6d-node%.2d.txt",laenge, messungen, node);//speichert Ableitung
	sprintf(dateinamezeit,"Messungen/Zeiten/zeiten-laenge-%.4d-m-%.6d-cores-%.2d-node%.2d.txt",laenge, messungen, anzahlcores, node);//speichert Ableitung
	bootstrapalledatei=fopen(dateinamebootstrapalle, "w+");
	ableitungdatei=fopen(dateinameableitung, "w");
	zeitdatei=fopen(dateinamezeit, "w");
	//Messen der zeit, die während des Programms vergeht, aus C-Kurs kopiert:
	struct timeval anfangbootstrap, endebootstrap, anfangprogramm, endeprogramm;
	double sec, usec, zeitgesbootstrap, summezeitgesbootstrap;
	summezeitgesbootstrap=0;//Zeit fuer alle Temperaturen insgesamt
	for (int n=0; n<temperaturzahl; n+=schritt){    //ueber alle gegebenen Temperaturen messen
		//printf("%d\n", n);
		sprintf(dateinamemessen,"Messungen/Messwerte/messung-laenge%.4d-m%.6d-t%.3d-node02.txt",laenge,messungen,n);//.2, damit alle dateinamengleich lang sind
		messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
		//gsl_rng_set(generator, seed);//initialisieren, bei jedem Durchlauf mit gleichem seed
		for(int core=0;core<anzahlcores;core+=1){
			generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
			gsl_rng_set(generatoren[core], seed+core);
		}
		gettimeofday(&anfangbootstrap, NULL);
		for(int len=0;len<12;len+=1){//Fuer verschiedene l blocking und bootstrapping durchfuehren
			l=blocklenarray[len];
			//printf("%d\t%d\n", n, l);
			if((blockarray=(double*)malloc(sizeof(double)*messungen/l))==NULL){//zum Speichern der Blocks, prüft, ob Speicherplatz richitg bereitgestellt wurde
				printf("Fehler beim Allokieren der Blocklaengen!\n");
				return (-1);
			};
			r=4*messungen;//Anzahl an Replikas, die beim Bootstrappen erzeugt werden
			blocks_generieren(l, messungen, 2, 6, blockarray, messdatei);//blocking
			//Vergleich bootstrapping mit und ohne parallelisierung
			bootstrap(l, r, messungen, temperaturarray[n], blockarray, generatoren,bootstrapalledatei);//bootstrapping
			//bootstrapohnepar(l, r, messungen, temperaturarray[n], blockarray, generatoren[0],bootstrapalledatei);//bootstrapping
			free(blockarray);
		}
		gettimeofday(&endebootstrap, NULL);
		sec= (double)(endebootstrap.tv_sec-anfangbootstrap.tv_sec);
		usec= (double)(endebootstrap.tv_usec-anfangbootstrap.tv_usec);
		zeitgesbootstrap=sec+1e-06*usec;
		summezeitgesbootstrap+=zeitgesbootstrap;
		fprintf(zeitdatei, "1.0\t%f\t%f\t%f\n", temperaturarray[n], (double)messungen, zeitgesbootstrap);
		printf("bei T=%f hat das Bootstrapping %f Sekunden gebraucht\n", temperaturarray[n], zeitgesbootstrap); 
		fclose(messdatei);
	}
	ableitung(128, temperaturzahl/schritt*12, 6, 3,4,5,1, bootstrapalledatei, ableitungdatei);
	gettimeofday(&endeprogramm, NULL);
	sec= (double)(endeprogramm.tv_sec-anfangprogramm.tv_sec);
	usec= (double)(endeprogramm.tv_usec-anfangprogramm.tv_usec);
	printf("Insgesamt hat das Botstrapping %f Sekunden gebraucht\n", summezeitgesbootstrap);
	fprintf(zeitdatei, "2.0\t-1.0\t%f\t%f\n",(double)messungen, summezeitgesbootstrap); 
	
	fclose(bootstrapalledatei);
	fclose(ableitungdatei);
	fclose(zeitdatei);
	free(temperaturarray);	
	//gsl_rng_free(generator);//free, close: zum Verhindern von Speicherproblemen
	for(int core=0;core<anzahlcores;core+=1){
		gsl_rng_free(generatoren[core]);
	}
	free(generatoren);
	return 0;
}
