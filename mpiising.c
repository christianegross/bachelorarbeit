//Christiane, start 10.07.20
//Misst Observablen Hamiltonian, Akzeptanz, Magnetisierung, und zweite und vierte Potenz der Magnetisierung bei verschiedenen Temperaturen
//Dazu erst thermalisierung des Gitters und dann Messungen
//Ueber Messergebnisse werden sowohl Standardschaetzer als auch Schaetzer mit Bootstrapping und Blocking von Mittelwert und Standardabweichung gebildet
//mit MPI parallelisiert


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <mpi.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "mpifunktionen.h"

int main(int argc, char **argv){
	MPI_Init(&argc, &argv);//Nur einmal im Pogramm moeglich, ganz am Anfang, damit alle Prozesse alle Variablen haben
	int myrank, anzahlprozesse;//fuer Aufteilungen benoetigt
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &anzahlprozesse);
	int laenge=anzahlprozesse*10;
	int schritt=10;
	if (argc>=3){//Auch Zuweisung ueber Kommandozeile moeglich
		//anzahlprozesse=atoi(argv[1]);
		laenge=atoi(argv[1]);
		schritt=atoi(argv[2]);
	}
	FILE *gitterdatei, *messdatei, *mitteldatei, *thermdatei;//Zum Speichern der Gitter, Messungen, naiven Mittelwerte, Messungen waehrend des thermalisierens
	char dateinamemittel[100], dateinametherm[100], dateinamemessen[100];//Dateinamen fuer Mittelwert, Gitter, Messungen
	double mittelwertmag, varianzmag, mittelwertakz, varianzakz;//fuer naive Fehler
	int messungen=10000;
	int N0;//Thermalisierungsschritte, je nach Temperatur unterschiedlich
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int temperaturzahl=650;//Temperaturen, bei denen gemessen wird
	int starttemp=0;
	int endtemp=temperaturzahl;
	double *temperaturarray;//Speichert verwendete Temperaturen
	if((temperaturarray=(double*)malloc(sizeof(double)*temperaturzahl))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren der Temperaturen!\n");
		return (-1);
	}
	for (int i=0; i<temperaturzahl;i++){//Temperaturarray intalisieren
		//~ //genaue Messung der Magnetisierung:
		//~ if((i<20)){temperaturarray[i]=0.05+i*0.1;}
		//~ if((i>=20)&&(i<48)){temperaturarray[i]=2.0+0.008*(i-20);}
		//~ if((i>=48)&&(i<136)){temperaturarray[i]=2.224+0.002*(i-48);}
		//~ if((i>=136)&&(i<181)){temperaturarray[i]=2.4+0.008*(i-136);}
		//~ if((i>=181)){temperaturarray[i]=2.76+0.032*(i-181);}
		//Um Hamiltonian vergleichen zu koennen
		if (i<50){temperaturarray[i]=i*0.02+0.02;}
		if ((i>=50)&&(i<300)){temperaturarray[i]=(i-50)*0.01+1;}
		if ((i>=300)&&(i<350)){temperaturarray[i]=(i-300)*0.02+3.5;}
		if ((i>=350)&&(i<450)){temperaturarray[i]=(i-350)*0.06+4.5;}
		if ((i>=450)&&(i<500)){temperaturarray[i]=(i-450)*0.2+10.01;}
		if ((i>=500)&&(i<550)){temperaturarray[i]=exp(2.996+(i-500)*0.032);}
		if ((i>=550)&&(i<600)){temperaturarray[i]=exp(4.605+(i-550)*0.046);}
		if ((i>=600)&&(i<650)){temperaturarray[i]=exp(6.907+(i-600)*0.046);}
	}
	gsl_rng **generatoren;//mehrere generatoren fuer parallelisierung
	if ((generatoren=(gsl_rng**)malloc(anzahlprozesse*sizeof(gsl_rng**)))==NULL){
		printf("Fehler beim Allokieren der Generatoren\n");
	}
	for(int core=0;core<anzahlprozesse;core+=1){//generatoren zuweisen und initialisieren
		generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(generatoren[core], seed+core);
	}
	sprintf(dateinamemittel, "Messungen/MPIMessungen/messenmittel-l%.4d-m-%.6d-proz%.2d-sch-%.2d.txt",laenge, messungen, anzahlprozesse, schritt);
	//In Dateien wird nur von einem Prozess geschrieben, daher nur ein Prozessmit "w" oeffnen, "w" kreiert Datei, deshalb warten, dass Datei sicher existiert, bis andere Prozesse Datei oeffnen
	if(myrank==0){
	mitteldatei=fopen(dateinamemittel, "w+");
	thermdatei=fopen("mpitestthermdummy.txt", "w+");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myrank!=0){
	mitteldatei=fopen(dateinamemittel, "r");
	thermdatei=fopen("mpitestthermdummy.txt", "r");
	}
	char gitter[laenge*laenge];//Initialisieren nur in einem Prozess, danach auf alle verteilen
	if(myrank==0){
		initialisierenhomogen(gitter, laenge);
	}
	MPI_Bcast(gitter, laenge*laenge, MPI_CHAR, 0, MPI_COMM_WORLD);
	for (int n=starttemp; n<endtemp; n+=schritt){    //ueber alle gegebenen Temperaturen messen
		if ((2<temperaturarray[n])&&(temperaturarray[n]<3)){N0=10000;}//In der Naehe des kritischen Punktes mehr Thermalisierungsschritte notwendig
		if ((2.25<temperaturarray[n])&&(temperaturarray[n]<2.4)){N0=20000;}
		if((2>=temperaturarray[n])||(3<=temperaturarray[n])) {N0=5000;}
		sprintf(dateinametherm,"Messungen/MPIMessungen/thermalisierung-laenge%.4d-m%.6d-t%.3d-proz%.2d.txt",laenge,messungen,n, anzahlprozesse);//.2, damit alle dateinamengleich lang sind
		sprintf(dateinamemessen,"Messungen/MPIMessungen/messung-laenge%.4d-m%.6d-t%.3d-proz%.2d.txt",laenge,messungen,n, anzahlprozesse);//.2, damit alle dateinamengleich lang sind
		//In Dateien wird nur von einem Prozess geschrieben, daher nur ein Prozessmit "w" oeffnen, "w" kreiert Datei, deshalb warten, dass Datei sicher existiert, bis andere Prozesse Datei oeffnen
		if(myrank==0){
			gitterdatei = fopen(dateinametherm, "w+");//Zum speichern der thermalisierten Gitter
			messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(myrank!=0){
			gitterdatei = fopen(dateinametherm, "r");
			messdatei = fopen(dateinamemessen, "r");
		}
		printf("%d\t%d\t", n, N0);
		//Erst Gitter thermalisieren, Messwerte davon nicht benoetigt
		thermalisierenmpi(N0, laenge, temperaturarray[n], j , gitter, thermdatei, gitterdatei, generatoren);
		//Danach Messungen, Ergebnisse davon verwenden
		messenmpi(messungen, laenge, temperaturarray[n], j, gitter, messdatei, generatoren);
		if(myrank==0){//Naive Auswertung, I/O benoetigt, daher nur ein Prozess
			mittelwertakz=mittelwertberechnungnaiv(messdatei, messungen, 1, 6);
			varianzakz=varianzberechnungnaiv(messdatei, messungen, mittelwertakz, 1, 6);
			mittelwertmag=mittelwertberechnungnaiv(messdatei, messungen, 2, 6);
			varianzmag=varianzberechnungnaiv(messdatei, messungen, mittelwertmag, 2, 6);
			fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", (double)laenge, temperaturarray[n],j/temperaturarray[n], mittelwertakz, varianzakz, mittelwertmag, varianzmag, temperaturarray[n]/j);
			printf("%f\t%f\n", mittelwertmag, varianzmag);
		}
		fclose(gitterdatei);
		fclose(messdatei);
	}
	
	//Bootstrapping der Ergebnisse, I/O benoetigt, daher nur ein Prozess
	if (myrank==0){
		int r;//replicas, die gezogen werden
		FILE *bootstrapalledateiakz, *bootstrapalledateimag, *bootstrapalledateimqu, *bootstrapalledateiham/*, *zeitdatei*/;//benoetigte Dateien zur Ausgabe
		char dateinamebootstrapalleakz[150], dateinamebootstrapallemag[150], dateinamebootstrapallemqu[150], dateinamebootstrapalleham[150]/*, dateinamezeit[150]*/;//Um Dateien mit Variablen benennen zu koennen
		int l;//Laenge der Blocks
		double *blockarray;//Zum Speichern der geblockten Messwerte
		double blocklenarray[10]={/*1,2,4,8,16,32, 64,*/128, 256, 384, 512, 640, 758, 876, 1024, 1280, 1536};//Blocklaengen, bei denen gemessen wird
		sprintf(dateinamebootstrapalleakz,"Messungen/Bootstrapges/bootstrapalle-akzeptanz-l%.4d-m-%.6d-proz%.2d-sch-%.2d.txt",laenge, messungen, anzahlprozesse, schritt);//speichert Mitteelwerte aus Bootstrap
		sprintf(dateinamebootstrapallemag,"Messungen/Bootstrapges/bootstrapalle-magnetisierung-l%.4d-m-%.6d-proz%.2d-sch-%.2d.txt",laenge, messungen, anzahlprozesse, schritt);//speichert Mitteelwerte aus Bootstrap
		sprintf(dateinamebootstrapallemqu,"Messungen/Bootstrapges/bootstrapalle-magquad-l%.4d-m-%.6d-proz%.2d-sch-%.2d.txt",laenge, messungen, anzahlprozesse, schritt);//speichert Mitteelwerte aus Bootstrap
		sprintf(dateinamebootstrapalleham,"Messungen/Bootstrapges/bootstrapalle-hamiltonian-l%.4d-m-%.6d-proz%.2d-sch-%.2d.txt",laenge, messungen, anzahlprozesse, schritt);//speichert Mitteelwerte aus Bootstrap
		bootstrapalledateiakz=fopen(dateinamebootstrapalleakz, "w+");
		bootstrapalledateimag=fopen(dateinamebootstrapallemag, "w+");
		bootstrapalledateimqu=fopen(dateinamebootstrapallemqu, "w+");
		bootstrapalledateiham=fopen(dateinamebootstrapalleham, "w+");
		
		for (int n=starttemp; n<endtemp; n+=schritt){    //ueber alle gegebenen Temperaturen Mittelwerte bilden
			sprintf(dateinamemessen,"Messungen/MPIMessungen/messung-laenge%.4d-m%.6d-t%.3d-proz%.2d.txt",laenge,messungen,n, anzahlprozesse);//.2, damit alle dateinamengleich lang sind
			messdatei = fopen(dateinamemessen, "r");//Zum Einlesen der Messdaten
			for(int len=0;len<10;len+=1){//Fuer verschiedene l blocking und bootstrapping durchfuehren
				l=blocklenarray[len];
				//printf("%d\t%d\n", n, l);
				if((blockarray=(double*)malloc(sizeof(double)*messungen/l))==NULL){//zum Speichern der Blocks, prüft, ob Speicherplatz richitg bereitgestellt wurde
					printf("Fehler beim Allokieren der Blocklaengen!\n");
					return (-1);
				};
				r=4*messungen;//Anzahl an Replikas, die beim Bootstrappen erzeugt werden
				//akzeptanzrate
				blocks_generieren(l, messungen, 1, 6, blockarray, messdatei);//blocking
				bootstrapohnepar(l, r, messungen, temperaturarray[n], blockarray, generatoren[0],bootstrapalledateiakz);//bootstrapping
				//magnetisierung
				blocks_generieren(l, messungen, 2, 6, blockarray, messdatei);//blocking
				bootstrapohnepar(l, r, messungen, temperaturarray[n], blockarray, generatoren[0],bootstrapalledateimag);//bootstrapping
				//magnetisierungquadrat
				blocks_generieren(l, messungen, 3, 6, blockarray, messdatei);//blocking
				bootstrapohnepar(l, r, messungen, temperaturarray[n], blockarray, generatoren[0],bootstrapalledateimqu);//bootstrapping
				//hamiltonian
				blocks_generieren(l, messungen, 5, 6, blockarray, messdatei);//blocking
				bootstrapohnepar(l, r, messungen, temperaturarray[n], blockarray, generatoren[0],bootstrapalledateiham);//bootstrapping
				free(blockarray);
			}
			fclose(messdatei);
		}
	}
	for(int core=0;core<anzahlprozesse;core+=1){
		gsl_rng_free(generatoren[core]);
	}
	free(generatoren);
	fclose(mitteldatei);
	printf("\n");
	
	MPI_Finalize();//benoetigt, damit Programm beendet werden kann
}

