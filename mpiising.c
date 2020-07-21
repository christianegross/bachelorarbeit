//Christiane, start 10.07.20
//Messungen mit MPI-Parallelisierung


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <mpi.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "mpifunktionen.h"

int main(int argc, char **argv){
	MPI_Init(&argc, &argv);
	//int anzahlprozesse;//=1;
	int myrank, anzahlprozesse;//, opened=1, *isopened;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &anzahlprozesse);
	int laenge=anzahlprozesse*10;
	int schritt=10;
	if (argc>=3){
		//anzahlprozesse=atoi(argv[1]);
		laenge=atoi(argv[1]);
		schritt=atoi(argv[2]);
	}
	FILE *gitterdatei, *messdatei, *mitteldatei, *thermdatei;
	char dateinamemittel[100], dateinametherm[100], dateinamemessen[100];
	double mittelwertmag, varianzmag, mittelwertakz, varianzakz;//fuer naive Fehler
	int messungen=10368;
	int N0;
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int temperaturzahl=210;//Temperaturen, beid enen gemessen wird
	int starttemp=0;
	int endtemp=temperaturzahl;
	double *temperaturarray;
	if((temperaturarray=(double*)malloc(sizeof(double)*temperaturzahl))==NULL){//speichert verwendete Temperaturen, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren der Temperaturen!\n");
		return (-1);
	}
	for (int i=0; i<temperaturzahl;i++){//Temperaturarray intalisieren
		//genaue Messung der Magnetisierung:
		if((i<20)){temperaturarray[i]=0.05+i*0.1;}
		if((i>=20)&&(i<48)){temperaturarray[i]=2.0+0.008*(i-20);}
		if((i>=48)&&(i<136)){temperaturarray[i]=2.224+0.002*(i-48);}
		if((i>=136)&&(i<181)){temperaturarray[i]=2.4+0.008*(i-136);}
		if((i>=181)){temperaturarray[i]=2.76+0.032*(i-181);}
	}
	gsl_rng **generatoren;//mehrere generatoren fuer parallelisierung
	if ((generatoren=(gsl_rng**)malloc(anzahlprozesse*sizeof(gsl_rng**)))==NULL){
		printf("Fehler beim Allokieren der Generatoren\n");
	}
	for(int core=0;core<anzahlprozesse;core+=1){
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
	int gitter[laenge*laenge];
	if(myrank==0){
		initialisierenhomogen(gitter, laenge);
	}
	MPI_Bcast(gitter, laenge*laenge, MPI_INT, 0, MPI_COMM_WORLD);
	for (int n=starttemp; n<endtemp; n+=schritt){    //ueber alle gegebenen Temperaturen messen
		if ((2<temperaturarray[n])&&(temperaturarray[n]<3)){N0=30000;}
		if ((2.25<temperaturarray[n])&&(temperaturarray[n]<2.4)){N0=100000;}
		if((2>=temperaturarray[n])||(temperaturarray[n]<=3)) {N0=5000;}
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
		printf("%d\t", n);
		thermalisierenmpi(N0, laenge, anzahlprozesse, myrank, temperaturarray[n], j , gitter, thermdatei, gitterdatei, generatoren);
		messenmpi(messungen, laenge, anzahlprozesse, myrank, temperaturarray[n], j, gitter, messdatei, generatoren);
		if(myrank==0){//Naive Auswertung
			mittelwertakz=mittelwertberechnungnaiv(messdatei, messungen, 1, 6);
			varianzakz=varianzberechnungnaiv(messdatei, messungen, mittelwertakz, 1, 6);
			mittelwertmag=mittelwertberechnungnaiv(messdatei, messungen, 2, 6);
			varianzmag=varianzberechnungnaiv(messdatei, messungen, mittelwertmag, 2, 6);
			fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", (double)laenge, temperaturarray[n],j/temperaturarray[n], mittelwertakz, varianzakz, mittelwertmag, varianzmag, temperaturarray[n]/j);
		}
		fclose(gitterdatei);
		fclose(messdatei);
	}
	
	//Bootstrapping der Ergebnisse
	if (myrank==0){
		int r;
		FILE *bootstrapalledateiakz, *bootstrapalledateimag, *bootstrapalledateimqu, *bootstrapalledateiham, *zeitdatei;//benoetigte Dateien zur Ausgabe
		char dateinamebootstrapalleakz[150], dateinamebootstrapallemag[150], dateinamebootstrapallemqu[150], dateinamebootstrapalleham[150], dateinamezeit[150];//Um Dateien mit Variablen benennen zu koennen
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
		
		for (int n=starttemp; n<endtemp; n+=schritt){    //ueber alle gegebenen Temperaturen messen
			sprintf(dateinamemessen,"Messungen/MPIMessungen/messung-laenge%.4d-m%.6d-t%.3d-proz%.2d.txt",laenge,messungen,n, anzahlprozesse);//.2, damit alle dateinamengleich lang sind
			messdatei = fopen(dateinamemessen, "r");//Zum Speichern der Messdaten
			for(int len=0;len<10;len+=1){//Fuer verschiedene l blocking und bootstrapping durchfuehren
				l=blocklenarray[len];
				//printf("%d\t%d\n", n, l);
				if((blockarray=(double*)malloc(sizeof(double)*messungen/l))==NULL){//zum Speichern der Blocks, prüft, ob Speicherplatz richitg bereitgestellt wurde
					printf("Fehler beim Allokieren der Blocklaengen!\n");
					return (-1);
				};
				r=4*messungen;//Anzahl an Replikas, die beim Bootstrappen erzeugt werden
				//akzaptanzrate
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
	
	MPI_Finalize();
}

