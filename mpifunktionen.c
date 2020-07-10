//Christiane, start 07.07.20
//Funktionen, die fuer die Parallelisierung mit MPI benoetigt werden.


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <mpi.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "auswertungsfunktionen.h"

int deltahmpi(int d1, int d2, int laenge, int teillaenge, int *untergitter, int *nachbarunten, int *nachbaroben){
	//berechnet Energieänderung bei Flip des Spins an position d1, d2
	//Dafuer jeweils Nachbarn noetig, nachbarn rechts links immer in Untergitter, aber oben/unten kann Teil eines anderen Gitters sein, deshalb Abgleich mit Nachbargittern
	int oben, unten, delta;
	if (d1==0){
		oben=nachbaroben[d2];
		unten=untergitter[laenge+d2];
	}
	if (d1==teillaenge-1){
		unten=nachbarunten[d2];
		oben=untergitter[(d1-1)*laenge+d2];
	}
	if ((d1!=teillaenge-1)&&(d1!=0)){
		unten=untergitter[(d1+1)*laenge+d2];
		oben=untergitter[(d1-1)*laenge+d2];
	}
	delta=2*untergitter[laenge*d1+d2]*(oben//oben
							  +unten//unten
							  +untergitter[laenge*d1+((d2-1+laenge)%(laenge))]//links
							  +untergitter[laenge*d1+((d2+1)%(laenge))]);//rechts
	return delta;
}

void initialisierenhomogen(int *gitter, int laenge){
	for (int i=1; i<laenge*laenge; i+=1){
		gitter[i]=1;
	}
}
	
//tryflip, einlesen, ausgabe hamiltonian: kopieren aus schon bekannten Funktionen, keine Aenderung noetig
void einlesen(int *gitter, int laenge, FILE *datei){
	//liest gitter, das in datei geschrieben wurde, wieder in gitter ein
	rewind(datei);
	int d1, d2, error;
	for (int n=0; n<laenge*laenge; n+=1){
		error= fscanf(datei, "%d \t %d \t %d \n", &d1, &d2, &gitter[n]);//Einlesen von Zeile, Spalte und Wert
		//Zeile und Spalte nicht benötigt, da Zuordnung desselbe wie beim Einlesen
		//Fehleranfällig? Doch über d1, d2 zuordnen?
		if (error<0){
			printf("Fehler beim Einlesen!\n");
			break;
		}
	}
}

void ausgabe(int *gitter, int laenge, FILE *datei){
	//gibt ein laenge *laenge quadratisches Gitter in datei aus
	for (int d1=0; d1<laenge; d1++){//geht in erster dimension durch
		for (int d2=0; d2<laenge; d2++){//geht in zweiter dimension durch
			fprintf(datei, "%3d\t%3d\t%d\n",d1, d2, gitter[laenge*d1+d2]);//Gibt Zeile, Spalte und Wert an
		}
	}
}

double hamiltonian(int *gitter, int laenge, double j){
	//berechnet den hamiltonian eines laenge*laenge quadratischen Gitters eines Ising Modells mit periodischen Randbedingungen
	double H=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			H+=(double)gitter[laenge*d1+d2]*(gitter[laenge*d1+((d2+1)%(laenge))]//Bond mit rechtem Nachbar
											+gitter[laenge*((d1+1)%(laenge))+d2]);//Bond mit unterem Nachbar
			//~ H+=(double)gitter[laenge*d1+d2]*gitter[laenge*((d1+1)%(laenge))+d2];//Bond mit unterem Nachbar
			//Von allen Feldern rechts und unten berücksichtig, periodische Randbedingungen->alles berücksichtigt
		}
	}
	return -j*H;	
}

extern inline int tryflip(gsl_rng *generator, double wahrscheinlichkeit){
	//versucht, den spin an position d1, d2 umzukehren nach Metropolis-Algorithmus
	double random=gsl_rng_uniform(generator);
	if (random<wahrscheinlichkeit){ 
	//if accepted return 1
		return 1;
	}
	else{
	//if rejected return 0
		return 0;
	}
	return -1;
}

double mittelwertberechnungnaiv(FILE *messdatei, int messungen, const int spalte, const int spalten){
	//berechnet den Mittelwert aus spalte aus den messungen in messdatei, messdatei muss nur aus doubles in spalten bestehen
	double summe=0;//Speichert Summe über Messungen
	double einwert=0;//Speichert einen ausgelesenen Wert
	double ergebnisarray[spalten];//speichert alle messungen
	rewind(messdatei);//sichergehen, dass alle Messdaten verwendet werden
	for (int messung=0; messung<messungen; messung+=1){//Mittelwert über Messung bilden
		for (int i=0; i<spalten; i+=1){
			fscanf(messdatei, "%le", &ergebnisarray[i]);//scannt einzelne doubles
			if (i==spalten-1){fscanf(messdatei, "\n");}//sorgt fuer Zeilenumbruch
		}
		einwert=ergebnisarray[spalte];//wählt korrekte messung aus
		summe+=einwert;
	}
	return summe/(double)messungen;
}

double varianzberechnungnaiv(FILE *messdatei, int messungen, double mittelwert, const int spalte, const int spalten){
	//berechnet die varianz über die gegebenen Messungen in messdatei mit mittelwert mittels Standardschaetzer, messdatei muss nur aus doubles in spalten bestehen
	double summe=0;//Speichert Summe über Messungen
	double einwert=0;//Speichert einen ausgelesenen Wert
	double ergebnisarray[spalten];//speichert alle messungen
	rewind(messdatei);//sichergehen, dass alle Messdaten verwendet werden
	for (int messung=0; messung<messungen; messung+=1){//Mittelwert über Messung bilden
		for (int i=0; i<spalten; i+=1){
			fscanf(messdatei, "%le", &ergebnisarray[i]);
			if (i==spalten-1){fscanf(messdatei, "\n");}
		}
		einwert=ergebnisarray[spalte];//wählt korrekte messung aus
		summe+=(einwert-mittelwert)*(einwert-mittelwert);
	}
	return sqrt(summe/((double)messungen-1));
}

double sweepmpi(int laenge, FILE *ausgabedatei, gsl_rng **generatoren, int anzproz, int myrank, int teillaenge, int *untergitter, int *nachbarunten, int *nachbaroben, double* wahrscheinlichkeiten, double j, double H){
	int delta;//potentielle Energieaenderung
	double ergebnisse[3]={0,0,0};//Nach sweep benutzen, um H, akz, mag aus reduzierung abzuspeichern
	double ergebnisselokal[3]={0,0,0};//Lokale Ergebnisse fuer H, akz, mag speichern
	//Untergitter durchgehen, an jedem Punkt Metropolis-Update
	for(int d1=0;d1<teillaenge;d1+=1){
		for(int d2=0;d2<laenge;d2+=1){
			delta=deltahmpi(d1, d2, laenge, teillaenge, untergitter, nachbarunten, nachbaroben);
			if (tryflip(generatoren[myrank], wahrscheinlichkeiten[delta/4+2])==1){
				//Wenn erfolgreich: Spinflip, H, akz aktualisieren
				untergitter[d1*laenge+d2]*=-1;
				ergebnisselokal[0]+=j*delta;//Zwischenvariable H
				ergebnisselokal[1]+=1;//Zwischenvariable akz
			}
			ergebnisselokal[2]+=(double)untergitter[d1*laenge+d2];//Berechnung mag, jeder Punkt noetig, nicht nur bei Spinflip
		}
	}
	//Gather nach jedem sweep noetig, damit nachbararrays gebildet werden können
	//ersetzen durch einzelne send/receive-Funktionen?
	MPI_Allreduce(ergebnisselokal, ergebnisse, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//H, akz, mag berechnen, indem Ergebnis aus jedem Teilgitter aufaddiert wird
	ergebnisse[0]+=H;//Bis jetzt wurde nur Veraenderung des Hamiltonians gemessen->jetzt auch aktualisierung auf Basis des vorherigen Wert
	if(myrank==0){
		double magnetisierung=(ergebnisse[2]>0?ergebnisse[2]/laenge/laenge:-ergebnisse[2]/laenge/laenge);//Betrag bilden
		double magquadrat=magnetisierung*magnetisierung;
		double magvier=magquadrat*magquadrat;
		fprintf(ausgabedatei, "%e\t%e\t%e\t%e\t%e\n",ergebnisse[1]/laenge/laenge, magnetisierung, magquadrat, magvier, ergebnisse[0]);//akz, mag, magquad, magvier, H
	}
	return ergebnisse[0];//Hamiltonian, fuer naechste Rechnung benoetit, zurueckgeben
			
}

void messenmpi(int messungen, int laenge, double T, double j, int *gitter, FILE *ausgabedatei, gsl_rng **generatoren){
	//int gitter[laenge*laenge];
	//for (int k=0;k<laenge*laenge;k+=1){gitter[k]=-1;};//Zum aisprobieren: Gitter am Anfang komplett homogen
	//einlesen(gitter, laenge, einlesedatei);
	double H=hamiltonian(gitter, laenge, j);//Anfangswert berechnen
	double wahrscheinlichkeiten[5]={1,1,1,exp(-4*j/T), exp(-8*j/T)};//Wahrscheinlichkeiten fuer Spinflip, muessen nur einmal berechnet werden
	if (j<0){
		wahrscheinlichkeiten[1]=wahrscheinlichkeiten[3];
		wahrscheinlichkeiten[0]=wahrscheinlichkeiten[4];
		wahrscheinlichkeiten[3]=1;
		wahrscheinlichkeiten[4]=1;
		}
	int anzproz, myrank;//Anzahl Prozessoren, meine Nummer
	MPI_Comm_size(MPI_COMM_WORLD, &anzproz);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int teillaenge=laenge/anzproz;
	int nachbaroben[laenge];
	int nachbarunten[laenge];
	int untergitter[laenge*teillaenge];		
	//eigenes Untergitter wird mit aktuellem Zustand des Gitters initialisiert. Nur einmal benoetigt, da danach Aenderungen von vorherigen sweep gespeichert sind.
	for (int i=0; i<teillaenge*laenge; i+=1){
		untergitter[i]=gitter[myrank*teillaenge*laenge+i];
	}
	for (int messung=0; messung<messungen; messung+=1){
		if(myrank==0){
			fprintf(ausgabedatei,"%e\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
			//printf("%e\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		}
		//Raender der benachbarten Untergitter werden gespeichert
		for (int i=0; i<laenge; i+=1){
			nachbaroben[i]=gitter[((teillaenge*myrank-1+laenge)%laenge)*laenge+i];
			nachbarunten[i]=gitter[((teillaenge*(myrank+1))%laenge)*laenge+i];
		}
		//eigenes Untergitter wird mit aktuellem Zustand des Gitters initialisiert.
		//for (int i=0; i<teillaenge*laenge; i+=1){
		//	untergitter[i]=gitter[myrank*teillaenge+i];
		//}
		H=sweepmpi(laenge,ausgabedatei, generatoren,  anzproz, myrank, teillaenge, untergitter, nachbarunten, nachbaroben,  wahrscheinlichkeiten,  j,  H);
		//Gitter aktualisieren nach sweep: In jedem Prozess Inhalte der einzelnen Untergitter in Gitter aneinanderreihen
		MPI_Allgather(untergitter, teillaenge, MPI_INT, gitter, teillaenge, MPI_INT, MPI_COMM_WORLD);
	}
		
}

void thermalisierenmpi(int messungen, int laenge, double T, double j, int *gitter, FILE *dummymessung, FILE *ausgabedatei, gsl_rng **generatoren){
	//int gitter[laenge*laenge];
	//for (int k=0;k<laenge*laenge;k+=1){gitter[k]=-1;};//Zum aisprobieren: Gitter am Anfang komplett homogen
	//einlesen(gitter, laenge, gitterdatei);
	double H=hamiltonian(gitter, laenge, j);//Anfangswert berechnen
	double wahrscheinlichkeiten[5]={1,1,1,exp(-4*j/T), exp(-8*j/T)};//Wahrscheinlichkeiten fuer Spinflip, muessen nur einmal berechnet werden
	if (j<0){
		wahrscheinlichkeiten[1]=wahrscheinlichkeiten[3];
		wahrscheinlichkeiten[0]=wahrscheinlichkeiten[4];
		wahrscheinlichkeiten[3]=1;
		wahrscheinlichkeiten[4]=1;
		}
	int anzproz, myrank;//Anzahl Prozessoren, meine Nummer
	MPI_Comm_size(MPI_COMM_WORLD, &anzproz);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int teillaenge=laenge/anzproz;
	int nachbaroben[laenge];
	int nachbarunten[laenge];
	int untergitter[laenge*teillaenge];		
	//eigenes Untergitter wird mit aktuellem Zustand des Gitters initialisiert. Nur einmal benoetigt, da danach Aenderungen von vorherigen sweep gespeichert sind.
	for (int i=0; i<teillaenge*laenge; i+=1){
		untergitter[i]=gitter[myrank*teillaenge*laenge+i];
	}
	for (int messung=0; messung<messungen; messung+=1){
		if(myrank==0){
			fprintf(dummymessung,"%e\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
			//printf("%e\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		}
		//Raender der benachbarten Untergitter werden gespeichert
		for (int i=0; i<laenge; i+=1){
			nachbaroben[i]=gitter[((teillaenge*myrank-1+laenge)%laenge)*laenge+i];
			nachbarunten[i]=gitter[((teillaenge*(myrank+1))%laenge)*laenge+i];
		}
		//eigenes Untergitter wird mit aktuellem Zustand des Gitters initialisiert.
		//for (int i=0; i<teillaenge*laenge; i+=1){
		//	untergitter[i]=gitter[myrank*teillaenge+i];
		//}
		H=sweepmpi(laenge,dummymessung, generatoren,  anzproz, myrank, teillaenge, untergitter, nachbarunten, nachbaroben,  wahrscheinlichkeiten,  j,  H);
		//Gitter aktualisieren nach sweep: In jedem Prozess Inhalte der einzelnen Untergitter in Gitter aneinanderreihen
		MPI_Allgather(untergitter, teillaenge, MPI_INT, gitter, teillaenge, MPI_INT, MPI_COMM_WORLD);
	}
	if(myrank==0){
		ausgabe(gitter, laenge, ausgabedatei);
	}
		
}





int main(int argc, char **argv){
	int anzahlprozesse=1;
	int laenge=anzahlprozesse*10;
	if (argc>=3){
		anzahlprozesse=atoi(argv[1]);
		laenge=atoi(argv[2]);
	}
	FILE*thermdatei=fopen("mpitestthermdummy.txt", "w+");
	FILE *gitterdatei, *messdatei, *mitteldatei;
	char dateinamemittel[100], dateinametherm[100], dateinamemessen[100];
	double mittelwertmag, varianzmag, mittelwertakz, varianzakz;//fuer naive Fehler
	int messungen=1000;
	int N0;
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int temperaturzahl=210;//Temperaturen, beid enen gemessen wird
	int starttemp=0;
	int endtemp=temperaturzahl;
	int schritt=5;
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
	gsl_rng **generatoren;//mehrere generaotren fuer parallelisierung
	if ((generatoren=(gsl_rng**)malloc(anzahlprozesse*sizeof(gsl_rng**)))==NULL){
		printf("Fehler beim Allokieren der Generatoren\n");
	}
	for(int core=0;core<anzahlprozesse;core+=1){
		generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(generatoren[core], seed+core);
	}
	sprintf(dateinamemittel, "Messungen/MPIMessungen/messenmittel-l%.4d-m-%.6d-proz%.2d-sch-%.2d.txt",laenge, messungen, anzahlprozesse, schritt);
	mitteldatei=fopen(dateinamemittel, "w+");
	int gitter[laenge*laenge];
	initialisierenhomogen(gitter, laenge);
	int myrank;//, opened=1, *isopened;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	for (int n=starttemp; n<endtemp; n+=schritt){    //ueber alle gegebenen Temperaturen messen
		if ((2<temperaturarray[n])&&(temperaturarray[n]<3)){N0=30000;}
		if ((2.25<temperaturarray[n])&&(temperaturarray[n]<2.4)){N0=100000;}
		if((2>=temperaturarray[n])||(temperaturarray[n]<=3)) {N0=5000;}
		sprintf(dateinametherm,"Messungen/MPIMessungen/thermalisierung-laenge%.4d-m%.6d-t%.3d-proz%.2d.txt",laenge,messungen,n, anzahlprozesse);//.2, damit alle dateinamengleich lang sind
		sprintf(dateinamemessen,"Messungen/MPIMessungen/messung-laenge%.4d-m%.6d-t%.3d-proz%.2d.txt",laenge,messungen,n, anzahlprozesse);//.2, damit alle dateinamengleich lang sind
		if(myrank==0){
			gitterdatei = fopen(dateinametherm, "w+");//Zum speichern der thermalisierten Gitter
			messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
			//isopened=&opened;
			//printf("%d\t", *isopened);
		}
		//MPI_Bcast(isopened, 1, MPI_INT, 0, MPI_COMM_WORLD);//Stellt sicher, dass Datei, die von nivht-schreibenden Prozessen geöffnter wird, auch existiert
		MPI_Barrier(MPI_COMM_WORLD);
		if(myrank!=0){
			gitterdatei = fopen(dateinametherm, "r");//Zum speichern der thermalisierten Gitter
			messdatei = fopen(dateinamemessen, "r");//Zum Speichern der Messdaten
		}
		printf("%d\t", n);
		thermalisierenmpi(N0, laenge, temperaturarray[n], j , gitter, thermdatei, gitterdatei, generatoren);
		messenmpi(messungen, laenge, temperaturarray[n], j, gitter, messdatei, generatoren);
		if(myrank==0){
			mittelwertakz=mittelwertberechnungnaiv(messdatei, messungen, 1, 6);
			varianzakz=varianzberechnungnaiv(messdatei, messungen, mittelwertakz, 1, 6);
			mittelwertmag=mittelwertberechnungnaiv(messdatei, messungen, 2, 6);
			varianzmag=varianzberechnungnaiv(messdatei, messungen, mittelwertmag, 2, 6);
			fprintf(mitteldatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", (double)laenge, temperaturarray[n],j/temperaturarray[n], mittelwertakz, varianzakz, mittelwertmag, varianzmag, temperaturarray[n]/j);
		}
		fclose(gitterdatei);
		fclose(messdatei);
	}
	MPI_Finalize();
	for(int core=0;core<anzahlprozesse;core+=1){
		gsl_rng_free(generatoren[core]);
	}
	free(generatoren);
	fclose(mitteldatei);
	printf("\n");
}


//christiane@christiane-VirtualBox2:17:07:41:~/Bachelorarbeit/bachelorarbeit$ mpicc -g -Wall -o mpifunktionen mpifunktionen.c -lgsl -lm -lgslcblas
//christiane@christiane-VirtualBox2:17:08:52:~/Bachelorarbeit/bachelorarbeit$ mpirun -n 10 ./mpifunktionen
