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
	int laenge=480;//laenge der verwendeten Gitter
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int N01=10000;//sweeps beim ersten Thermalisieren
	int N0=5000;//benoetigte sweeps zum Thermalisieren
	int messungen=10000;//pro temperatur, zweierpotenz um blocken einfacher zu machen
	int r;//Anzahl an samples für den Bootstrap
	FILE *gitterthermdatei, *messdatei, *mittelwertdatei, *dummydatei, *bootstrapalledateiakz, *bootstrapalledateimag, *bootstrapalledateiham, *ableitungdatei, *zeitdatei;//benoetigte Dateien zur Ausgabe
	int temperaturzahl=100;//Temperaturen, beid enen gemessen wird
	int schritt=1;//Wie viele Punkte werden gemessen?
	int node=2;//nodes auf vm, qbig
	char dateinametherm[150], dateinamemessen[150], dateinamemittel[150], dateinamebootstrapalleakz[150], dateinamebootstrapallemag[150], dateinamebootstrapalleham[150], dateinameableitung[150], dateinamezeit[150];//Um Dateien mit Variablen benennen zu koennen
	double mittelwertmag, varianzmag, mittelwertakz, varianzakz;//fuer naive Fehler
	double U, magquad, varmagquad, magvier, varmagvier;
	double *temperaturarray;
	if((temperaturarray=(double*)malloc(sizeof(double)*temperaturzahl))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren der Temperaturen!\n");
		return (-1);
	}
	for (int i=0; i<temperaturzahl;i++){//Temperaturarray intalisieren
		//MAgnetisierung zur genaueren Bestimmung des kritishen Punktes
		temperaturarray[i]=2.1+0.0025;
		//Für Akzeptanzrate, um bis 10.000 zu kommen, gemessen mit Schritt 2
		//~ if (i<50){temperaturarray[i]=i*0.02+0.02;}
		//~ if ((i>=50)&&(i<300)){temperaturarray[i]=(i-50)*0.01+1;}
		//~ if ((i>=300)&&(i<350)){temperaturarray[i]=(i-300)*0.02+3.5;}
		//~ if ((i>=350)&&(i<450)){temperaturarray[i]=(i-350)*0.06+4.5;}
		//~ if ((i>=450)&&(i<500)){temperaturarray[i]=(i-450)*0.2+10.01;}
		//~ if ((i>=500)&&(i<550)){temperaturarray[i]=exp(2.996+(i-500)*0.032);}
		//~ if ((i>=550)&&(i<600)){temperaturarray[i]=exp(4.605+(i-550)*0.046);}
		//~ if ((i>=600)&&(i<650)){temperaturarray[i]=exp(6.907+(i-600)*0.046);}
		//printf("%d\t%e\n", i, temperaturarray[i]);
	}
	printf("Ende array\n");
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
	
	sprintf(dateinamemittel,"Messungen/Mittelwerte/messenmittel-l%.4d-m-%.6d-node%.2d.txt",laenge, messungen, node);//speichert naive Mittelwerte
	sprintf(dateinamebootstrapalleakz,"Messungen/Bootstrapges/bootstrapalle-akzeptanz-l%.4d-m-%.6d-node%.2d.txt",laenge, messungen, node);//speichert Mitteelwerte aus Bootstrap
	sprintf(dateinamebootstrapallemag,"Messungen/Bootstrapges/bootstrapalle-magnetisierung-l%.4d-m-%.6d-node%.2d.txt",laenge, messungen, node);//speichert Mitteelwerte aus Bootstrap
	sprintf(dateinamebootstrapalleham,"Messungen/Bootstrapges/bootstrapalle-hamiltonian-l%.4d-m-%.6d-node%.2d.txt",laenge, messungen, node);//speichert Mitteelwerte aus Bootstrap
	sprintf(dateinameableitung,"Messungen/ableitung-akzeptanz-laenge-%.4d-m-%.6d-node%.2d.txt",laenge, messungen, node);//speichert Ableitung
	sprintf(dateinamezeit,"Messungen/Zeiten/zeiten-laenge-%.4d-m-%.6d-cores-%.2d-node%.2d.txt",laenge, messungen, anzahlcores, node);//speichert Ableitung
	mittelwertdatei=fopen(dateinamemittel, "w+");
	bootstrapalledateiakz=fopen(dateinamebootstrapalleakz, "w+");
	bootstrapalledateimag=fopen(dateinamebootstrapallemag, "w+");
	bootstrapalledateiham=fopen(dateinamebootstrapalleham, "w+");
	ableitungdatei=fopen(dateinameableitung, "w");
	zeitdatei=fopen(dateinamezeit, "w");
	//Messen der zeit, die während des Programms vergeht, aus C-Kurs kopiert:
	struct timeval anfangmessen, endemessen, anfangbootstrap, endebootstrap, anfangprogramm, endeprogramm;
	double sec, usec, zeitgesmessen, summezeitgesmessen, zeitgesbootstrap, summezeitgesbootstrap, zeitgesprogramm;
	summezeitgesmessen=0;//Zeit fuer alle Temperaturen insgesamt
	summezeitgesbootstrap=0;//Zeit fuer alle Temperaturen insgesamt
	gettimeofday(&anfangprogramm, NULL);
	//Thermalisierung des ersten Gitters, nicht ueber letztes verwendetes Gitter moeglich
	char gitter[laenge*laenge];
	initialisierung(gitter, laenge, seed);
	dummydatei=fopen("dummytherm.txt", "w");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird
	//gsl_rng_set(generator, seed);
	thermalisierenmehreregeneratoren(laenge, temperaturarray[0], j, seed, N01, gitter, dummydatei, generatoren);//Erstes Thermalisierens, Anzahl je nach Länge groesser machen
	fclose(dummydatei);
	for (int n=0; n<temperaturzahl; n+=schritt){    //ueber alle gegebenen Temperaturen messen
		printf("%d\n", n);
		sprintf(dateinametherm,"Messungen/ThermalisierteGitter/thermalisierung-laenge%.4d-m%.6d-t%.3d-node%.2d.txt",laenge,messungen,n, node);//.2, damit alle dateinamengleich lang sind
		sprintf(dateinamemessen,"Messungen/Messwerte/messung-laenge%.4d-m%.6d-t%.3d-node%.2d.txt",laenge,messungen,n, node);//.2, damit alle dateinamengleich lang sind
		gitterthermdatei = fopen(dateinametherm, "w+");//Zum speichern der thermalisierten Gitter
		messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
		//gsl_rng_set(generator, seed);//initialisieren, bei jedem Durchlauf mit gleichem seed
		for(int core=0;core<anzahlcores;core+=1){
			generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
			gsl_rng_set(generatoren[core], seed+core);
		}
		
		thermalisierenmehreregeneratoren(laenge, temperaturarray[n], j, seed, N0, gitter, gitterthermdatei, generatoren);
		gettimeofday(&anfangmessen, NULL);
		messenmehreregeneratoren(laenge, temperaturarray[n], j, messungen, gitter/*thermdatei*/, messdatei, generatoren);
		gettimeofday(&endemessen, NULL);
		sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
		usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
		zeitgesmessen=sec+1e-06*usec;
		summezeitgesmessen+=zeitgesmessen;
		//printf("bei T=%f haben %d Messungen %f Sekunden gebraucht\n", temperaturarray[n], messungen, zeitgesmessen);
		fprintf(zeitdatei, "0.0\t%f\t%f\t%f\n", temperaturarray[n], (double)messungen, zeitgesmessen);
		//Berechnung der naiven Standardfehler
		mittelwertakz=mittelwertberechnungnaiv(messdatei, messungen, 1, 6);
		varianzakz=varianzberechnungnaiv(messdatei, messungen, mittelwertakz, 1, 6);
		mittelwertmag=mittelwertberechnungnaiv(messdatei, messungen, 2, 6);
		varianzmag=varianzberechnungnaiv(messdatei, messungen, mittelwertmag, 2, 6);
		magquad=mittelwertberechnungnaiv(messdatei, messungen, 3, 6);
		magvier=mittelwertberechnungnaiv(messdatei, messungen, 4, 6);
		U=1-(magvier/3/magquad/magquad);
		fprintf(mittelwertdatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", (double)laenge, temperaturarray[n],j/temperaturarray[n], mittelwertakz, varianzakz, mittelwertmag, varianzmag, temperaturarray[n]/j, U);
		gettimeofday(&anfangbootstrap, NULL);
		for(int len=0;len<12;len+=1){//Fuer verschiedene l blocking und bootstrapping durchfuehren
			l=blocklenarray[len];
			//printf("%d\t%d\n", n, l);
			if((blockarray=(double*)malloc(sizeof(double)*messungen/l))==NULL){//zum Speichern der Blocks, prüft, ob Speicherplatz richitg bereitgestellt wurde
				printf("Fehler beim Allokieren der Blocklaengen!\n");
				return (-1);
			};
			r=4*messungen;//Anzahl an Replikas, die beim Bootstrappen erzeugt werden
			//~ //akzaptanzrate
			//~ blocks_generieren(l, messungen, 1, 6, blockarray, messdatei);//blocking
			//~ //Vergleich bootstrapping mit und ohne parallelisierung
			//~ bootstrap(l, r, messungen, temperaturarray[n], blockarray, generatoren,bootstrapalledateiakz);//bootstrapping
			//magnetisierung
			blocks_generieren(l, messungen, 2, 6, blockarray, messdatei);//blocking
			//Vergleich bootstrapping mit und ohne parallelisierung
			bootstrap(l, r, messungen, temperaturarray[n], blockarray, generatoren,bootstrapalledateimag);//bootstrapping
			//bootstrapohnepar(l, r, messungen, temperaturarray[n], blockarray, generator,bootstrapalledatei);//bootstrapping
			//~ //hamiltonian
			//~ blocks_generieren(l, messungen, 5, 6, blockarray, messdatei);//blocking
			//~ //Vergleich bootstrapping mit und ohne parallelisierung
			//~ bootstrap(l, r, messungen, temperaturarray[n], blockarray, generatoren,bootstrapalledateiham);//bootstrapping
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
	ableitung(128, temperaturzahl/schritt*12, 6, 3,4,5,1, bootstrapalledateimag, ableitungdatei);
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
	fclose(bootstrapalledateiakz);
	fclose(bootstrapalledateimag);
	fclose(bootstrapalledateiham);
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
