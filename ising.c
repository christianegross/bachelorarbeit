//Christiane, start 15.04.20
//Erste Implementierung Ising-Modell für Bachelorarbeit
//Versuch Metropolis-Algorithmus


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "messfunktionen.h"
#include "auswertungsfunktionen.h"

int main(int argc, char **argv){
	//benoetigte Variablen initialisieren
	int anzahlcores=1;
	omp_set_num_threads(anzahlcores);//Setzt die nummer an Kernen, die in den parallelen Regionen verwendet werden.
	int laenge=50;//laenge der verwendeten Gitter
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int N01=1000;//sweeps beim ersten Thermalisieren
	int N0=1000;//benoetigte sweeps zum Thermalisieren
	int messungen=10000;//pro temperatur, zweierpotenz um blocken einfacher zu machen
	int r;//Anzahl an samples für den Bootstrap
	FILE *gitterthermdatei, *messdatei, *mittelwertdatei, *dummydatei, *bootstrapalledatei, *ableitungdatei, *zeitdatei;//benoetigte Dateien zur Ausgabe
	int temperaturzahl=300;//Temperaturen, beid enen gemessen wird
	char dateinametherm[150], dateinamemessen[150], dateinamemittel[150], dateinamebootstrapalle[150], dateinameableitung[150], dateinamezeit[150];//Um Dateien mit Variablen benennen zu koennen
	double mittelwertmag, varianzmag, mittelwertakz, varianzakz;//fuer naive Fehler
	double *temperaturarray;
	if((temperaturarray=(double*)malloc(sizeof(double)*temperaturzahl))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren der Temperaturen!\n");
		return (-1);
	}
	for (int i=0; i<temperaturzahl;i++){//Temperaturarray intalisieren
		temperaturarray[i]=0.015*i+0.015;
	}
	int l;//Laenge der Blocks
	double *blockarray;//Zum Speichern der geblockten Messwerte
	double blocklenarray[12]={32, 64,128, 256, 384, 512, 640, 758, 876, 1024, 1280, 1536};//Blocklaengen, bei denen gemessen wird
	gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	sprintf(dateinamemittel,"Messungen/Mittelwerte/messenmittel-l%.4d-m-%.6d.txt",laenge, messungen);//speichert naive Mittelwerte
	sprintf(dateinamebootstrapalle,"Messungen/Bootstrapges/bootstrapalle-l%.4d-m-%.6d.txt",laenge, messungen);//speichert Mitteelwerte aus Bootstrap
	sprintf(dateinameableitung,"Messungen/ableitung-laenge-%.4d-m-%.6d.txt",laenge, messungen);//speichert Ableitung
	sprintf(dateinamezeit,"Messungen/zeiten-laenge-%.4d-m-%.6d-cores-%.2d.txt",laenge, messungen, anzahlcores);//speichert Ableitung
	mittelwertdatei=fopen(dateinamemittel, "w+");
	bootstrapalledatei=fopen(dateinamebootstrapalle, "w+");
	ableitungdatei=fopen(dateinameableitung, "w");
	zeitdatei=fopen(dateinamezeit, "w");
	//Messen der zeit, die während des Programms vergeht, aus C-Kurs kopiert:
	struct timeval anfangmessen, endemessen, anfangbootstrap, endebootstrap, anfangprogramm, endeprogramm;
	double sec, usec, zeitgesmessen, summezeitgesmessen, zeitgesbootstrap, summezeitgesbootstrap, zeitgesprogramm;
	summezeitgesmessen=0;//Zeit fuer alle Temperaturen insgesamt
	summezeitgesbootstrap=0;//Zeit fuer alle Temperaturen insgesamt
	gettimeofday(&anfangprogramm, NULL);
	//Thermalisierung des ersten Gitters, nicht ueber letztes verwendetes Gitter moeglich
	int gitter[laenge*laenge];
	initialisierung(gitter, laenge, seed);
	dummydatei=fopen("dummytherm.txt", "w");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird
	gsl_rng_set(generator, seed);
	thermalisieren(laenge, temperaturarray[0], j, seed, N01, gitter, dummydatei, generator);//Erstes Thermalisierens, Anzahl je nach Länge groesser machen
	fclose(dummydatei);
	for (int n=0; n<temperaturzahl; n+=1){    //ueber alle gegebenen Temperaturen messen
		printf("%d\n", n);
		sprintf(dateinametherm,"Messungen/ThermalisierteGitter/thermalisierung-laenge%.4d-m%.6d-t%.3d.txt",laenge,messungen,n);//.2, damit alle dateinamengleich lang sind
		sprintf(dateinamemessen,"Messungen/Messwerte/messung-laenge%.4d-m%.6d-t%.3d.txt",laenge,messungen,n);//.2, damit alle dateinamengleich lang sind
		gitterthermdatei = fopen(dateinametherm, "w+");//Zum speichern der thermalisierten Gitter
		messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
		gsl_rng_set(generator, seed);//initialisieren, bei jedem Durchlauf mit gleichem seed
		thermalisieren(laenge, temperaturarray[n], j, seed, N0, gitter, gitterthermdatei, generator);
		gettimeofday(&anfangmessen, NULL);
		messen(laenge, temperaturarray[n], j, messungen, gitterthermdatei, messdatei, generator);
		gettimeofday(&endemessen, NULL);
		sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
		usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
		zeitgesmessen=sec+1e-06*usec;
		summezeitgesmessen+=zeitgesmessen;
		//printf("bei T=%f haben %d Messungen %f Sekunden gebraucht\n", temperaturarray[n], messungen, zeitgesmessen);
		fprintf(zeitdatei, "0.0\t%f\t%f\t%f\n", temperaturarray[n], (double)messungen, zeitgesmessen);
		mittelwertakz=mittelwertberechnungnaiv(messdatei, messungen, 1, 3);
		varianzakz=varianzberechnungnaiv(messdatei, messungen, mittelwertakz, 1, 3);
		mittelwertmag=mittelwertberechnungnaiv(messdatei, messungen, 2, 3);
		varianzmag=varianzberechnungnaiv(messdatei, messungen, mittelwertmag, 2, 3);
		fprintf(mittelwertdatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", (double)laenge, temperaturarray[n],j/temperaturarray[n], mittelwertakz, varianzakz, mittelwertmag, varianzmag);
		gettimeofday(&anfangbootstrap, NULL);
		for(int len=0;len<12;len+=1){//Fuer verschiedene l blocking und bootstrapping durchfuehren
			l=blocklenarray[len];
			//printf("%d\t%d\n", n, l);
			if((blockarray=(double*)malloc(sizeof(double)*messungen/l))==NULL){//zum Speichern der Blocks, prüft, ob Speicherplatz richitg bereitgestellt wurde
				printf("Fehler beim Allokieren der Blocklaengen!\n");
				return (-1);
			};
			r=1;//r=4*messungen;//Anzahl an Replikas, die beim Bootstrappen erzeugt werden
			blocks_generieren(l, messungen, 2, 3, blockarray, messdatei);//blocking
			//Vergleich bootstrapping mit und ohne parallelisierung
			bootstrap(l, r, messungen, temperaturarray[n], blockarray, generator,bootstrapalledatei);//bootstrapping
			//bootstrapohnepar(l, r, messungen, temperaturarray[n], blockarray, generator,bootstrapalledatei);//bootstrapping
			free(blockarray);
		}
		gettimeofday(&endebootstrap, NULL);
		sec= (double)(endebootstrap.tv_sec-anfangbootstrap.tv_sec);
		usec= (double)(endebootstrap.tv_usec-anfangbootstrap.tv_usec);
		zeitgesbootstrap=sec+1e-06*usec;
		summezeitgesbootstrap+=zeitgesbootstrap;
		fprintf(zeitdatei, "1.0\t%f\t%f\t%f\n", temperaturarray[n], (double)messungen, zeitgesbootstrap);
		//printf("bei T=%f hat das Bootstrapping %f Sekunden gebraucht\n", temperaturarray[n], zeitgesbootstrap); 
		fclose(messdatei);
		fclose(gitterthermdatei);
	}
	//ableitung(128, 30*12, 6, 3,4,5,1, bootstrapalledatei, ableitungdatei);
	gettimeofday(&endeprogramm, NULL);
	sec= (double)(endeprogramm.tv_sec-anfangprogramm.tv_sec);
	usec= (double)(endeprogramm.tv_usec-anfangprogramm.tv_usec);
	zeitgesprogramm=sec+1e-06*usec;
	printf("Insgesamt hat das Messen %f Sekunden gebraucht\n", summezeitgesmessen);
	printf("Insgesamt hat das Botstrapping %f Sekunden gebraucht\n", summezeitgesbootstrap);
	printf("Insgesamt hat das Programm %f Sekunden gebraucht\n", zeitgesprogramm);
	fprintf(zeitdatei, "2.0\t-1.0\t%f\t%f\n",(double)messungen, summezeitgesbootstrap); 
	fprintf(zeitdatei, "3.0\t-1.0\t%f\t%f\n",(double)messungen, summezeitgesmessen); 
	fprintf(zeitdatei, "4.0\t-1.0\t%f\t%f\n",(double)messungen, zeitgesprogramm); 
	fclose(mittelwertdatei);
	fclose(bootstrapalledatei);
	fclose(ableitungdatei);
	fclose(zeitdatei);
	free(temperaturarray);
	gsl_rng_free(generator);//free, close: zum Verhindern von Speicherproblemen

	return 0;
}
