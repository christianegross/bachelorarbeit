//Christiane, start 07.07.20
//Funktionen, die fuer die Parallelisierung mit MPI benoetigt werden.


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <mpi.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen

int deltah(int d1, int d2, int laenge, int teillaenge, int *untergitter, int *nachbarunten, int *nachbaroben){
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
	
//tryflip, einlesen, hamiltonian: kopieren aus schon bekannten Funktionen, keine Aenderung noetig
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

double sweepmpi(int laenge, FILE *ausgabedatei, gsl_rng **generatoren, int anzproz, int myrank, int teillaenge, int *gitter, double* wahrscheinlichkeiten, double j, double H){
	int nachbaroben[laenge];
	int nachbarunten[laenge];
	int delta;//potentielle Energieaenderung
	int gitterelemente=0;//Ueberpruefen, ob alle Punkte behandelt werden
	//Untergitter schon als Parameter uebergeben, Verteilung durch scatter vor sweep/SChleife?
	int untergitter[laenge*teillaenge];
	//eigenes Untergitter wird mit aktuellem Zustand des Gitters initialisiert.
	for (int i=0; i<laenge; i+=1){
		nachbaroben[i]=gitter[((teillaenge*myrank-1+laenge)%laenge)*laenge+i];
		nachbarunten[i]=gitter[((teillaenge*myrank+1)%laenge)*laenge+i];
	}
	//Raender der benachbarten Untergitter werden gespeichert
	for (int i=0; i<teillaenge*laenge; i+=1){
		untergitter[i]=gitter[myrank*teillaenge+i];
	}
	double ergebnisse[3]={0,0,0};//Nach sweep benutzen, um H, akz, mag aus reduzierung abzuspeichern
	double ergebnisselokal[3]={0,0,0};//Lokale Ergebnisse fuer H, akz, mag speichern
	//Untergitter durchgehen, an jedem Punkt Metropolis-Update
	for(int d1=0;d1<teillaenge;d1+=1){
		for(int d2=0;d2<laenge;d2+=1){
			delta=deltah(d1, d2, laenge, teillaenge, untergitter, nachbarunten, nachbaroben);
			if (tryflip(generatoren[myrank], wahrscheinlichkeiten[delta/4+2])==1){
				//Wenn erfolgreich: Spinflip, H, akz aktualisieren
				untergitter[d1*laenge+d2]*=-1;
				ergebnisselokal[0]+=j*delta;//Zwischenvariable H
				ergebnisselokal[1]+=1;//Zwischenvariable akz
			}
			ergebnisselokal[2]+=(double)untergitter[d1*laenge+d2];//Berechnung mag, jeder Punkt noetig, nicht nur bei Spinflip
			gitterelemente+=1;//Ueberpruefung
		}
	}
	//Gather nach jedem sweep noetig, damit nachbararrays gebildet werden können
	//ersetzen durch einzelne send/receive-Funktionen?
	//Gitter aktualisieren nach sweep: In jedem Prozess Inhalte der einzelnen Untergitter in Gitter aneinanderreihen
	MPI_Allgather(untergitter, teillaenge, MPI_INT, gitter, teillaenge, MPI_INT, MPI_COMM_WORLD);
	MPI_Allreduce(ergebnisselokal, ergebnisse, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//H, akz, mag berechnen, indem Ergebnis aus jedem Teilgitter aufaddiert wird
	ergebnisse[0]+=H;//Bis jetzt wurde nur Veraenderung des Hamiltonians gemessen->jetzt auch aktualisierung auf Basis des vorherigen Wert
	//printf("%f\t", ergebnisse[2]/teillaenge/laenge);
	if(myrank==0){
	double magnetisierung=(ergebnisse[2]>0?ergebnisse[2]/laenge/laenge:-ergebnisse[2]/laenge/laenge);//Betrag bilden
	double magquadrat=magnetisierung*magnetisierung;
	double magvier=magquadrat*magquadrat;
	fprintf(ausgabedatei, "%e\t%e\t%e\t%e\t%e\n",ergebnisse[1]/laenge/laenge, magnetisierung, magquadrat, magvier, ergebnisse[0]);//akz, mag, magquad, magvier, H
	//printf("%d\t", gitterelemente);
	}
	return ergebnisse[0];//Hamiltonian, fuer naechste Rechnung benoetit, zurueckgeben
			
}

void messenmpi(int messungen, int laenge, double T, double j,/* FILE *einlesedatei, */FILE *ausgabedatei, gsl_rng **generatoren){
	int gitter[laenge*laenge];
	for (int k=0;k<laenge*laenge;k+=1){gitter[k]=-1;};//Zum aisprobieren: Gitter am Anfang komplett homogen
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
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &anzproz);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	int teillaenge=laenge/anzproz;
	for (int messung=0; messung<messungen; messung+=1){
		if(myrank==0){
			fprintf(ausgabedatei,"%e\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
			//printf("%e\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		}
		H=sweepmpi(laenge,ausgabedatei, generatoren,  anzproz, myrank, teillaenge, gitter,  wahrscheinlichkeiten,  j,  H);
	}
	MPI_Finalize();
		
}





int main(int argc, char **argv){
	FILE*datei=fopen("mpitest.txt", "w+");
	gsl_rng **generatoren;//mehrere generaotren fuer parallelisierung
	if ((generatoren=(gsl_rng**)malloc(10*sizeof(gsl_rng**)))==NULL){
		printf("Fehler beim Allokieren der Generatoren\n");
	}
	for(int core=0;core<10;core+=1){
		generatoren[core]=gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(generatoren[core], 5+core);
	}
	messenmpi(100, 50, 4.5, 1, datei, generatoren);
	for(int core=0;core<10;core+=1){
		gsl_rng_free(generatoren[core]);
	}
	free(generatoren);
	fclose(datei);
}

//christiane@christiane-VirtualBox2:17:07:41:~/Bachelorarbeit/bachelorarbeit$ mpicc -g -Wall -o mpifunktionen mpifunktionen.c -lgsl -lm -lgslcblas
//christiane@christiane-VirtualBox2:17:08:52:~/Bachelorarbeit/bachelorarbeit$ mpirun -n 10 ./mpifunktionen
