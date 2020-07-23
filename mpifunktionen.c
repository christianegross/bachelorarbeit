//Christiane, start 07.07.20
//Funktionen, die fuer die Parallelisierung mit MPI benoetigt werden.


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <mpi.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "mpifunktionen.h"

int deltahmpi(int d1, int d2, int laenge, int teillaenge, char *untergitter, char *nachbarunten, char *nachbaroben){
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

void initialisierenhomogen(char *gitter, int laenge){
	for (int i=0; i<laenge*laenge; i+=1){
		gitter[i]=1;
	}
}
	
//tryflip, einlesen, ausgabe hamiltonian: kopieren aus schon bekannten Funktionen, keine Aenderung noetig
void einlesen(char *gitter, int laenge, FILE *datei){
	//liest gitter, das in datei geschrieben wurde, wieder in gitter ein
	rewind(datei);
	int d1, d2, error;
	for (int n=0; n<laenge*laenge; n+=1){
		error= fscanf(datei, "%d \t %d \t %c \n", &d1, &d2, &gitter[n]);//Einlesen von Zeile, Spalte und Wert
		//Zeile und Spalte nicht benötigt, da Zuordnung desselbe wie beim Einlesen
		//Fehleranfällig? Doch über d1, d2 zuordnen?
		if (error<0){
			printf("Fehler beim Einlesen!\n");
			break;
		}
	}
}

void ausgabe(char *gitter, int laenge, FILE *datei){
	//gibt ein laenge *laenge quadratisches Gitter in datei aus
	for (int d1=0; d1<laenge; d1++){//geht in erster dimension durch
		for (int d2=0; d2<laenge; d2++){//geht in zweiter dimension durch
			fprintf(datei, "%3d\t%3d\t%c\n",d1, d2, gitter[laenge*d1+d2]);//Gibt Zeile, Spalte und Wert an
		}
	}
}

double hamiltonian(char *gitter, int laenge, double j){
	//berechnet den hamiltonian eines laenge*laenge quadratischen Gitters eines Ising Modells mit periodischen Randbedingungen
	double H=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			H+=(double)gitter[laenge*d1+d2]*(gitter[laenge*d1+((d2+1)%(laenge))]//Bond mit rechtem Nachbar
											+gitter[laenge*((d1+1)%(laenge))+d2]);//Bond mit unterem Nachbar
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

double sweepmpi(int laenge, FILE *ausgabedatei, gsl_rng **generatoren, int anzproz, int myrank, int teillaenge, char *untergitter, char *nachbarunten, char *nachbaroben, double* wahrscheinlichkeiten, double j, double H){
	//Fuehrt ein Metropolis-Update an jedem Punkt des gegeben Untergitters aus, durch Kombination aus mehreren Prozessen ein Metropolis-Update fuer jeden Punkt
	
	//printf("anfang sf in %d\n", myrank);
	int delta;//potentielle Energieaenderung
	double ergebnisse[3]={0,0,0};//Nach sweep benutzen, um H, akz, mag aus reduzierung abzuspeichern
	double ergebnisselokal[3]={0,0,0};//Lokale Ergebnisse fuer H, akz, mag speichern
	int sendbufferoben[laenge], sendbufferunten[laenge];
	//Untergitter durchgehen, an jedem Punkt Metropolis-Update
	//erst schwarze Punkte:
	
	//printf("anfang sweep schwarz in %d\n", myrank);
	for(int d1=0;d1<teillaenge;d1+=1){
		for(int d2=(teillaenge*anzproz+d1)%2;d2<laenge;d2+=2){
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
	
	//printf("sweep schwarz in %d\n", myrank);
	//Austausch Raender
	for(int i=0;i<laenge;i+=1){
		//rechter Rand
		sendbufferunten[i]=untergitter[(laenge*(teillaenge-1))+i];
		//linker Rand
		sendbufferoben[i]=untergitter[i];
	}
	//Erst Neue Raender nach oben senden/von unten empfangen
	MPI_Sendrecv(sendbufferoben, laenge, MPI_CHAR, (myrank-1+anzproz)%anzproz, 0, nachbarunten, laenge, MPI_CHAR, (myrank+1)%anzproz, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//Dann nach unten senden/von oben empfangen -> beide Male geschlossener Kreis
	MPI_Sendrecv(sendbufferunten, laenge, MPI_CHAR, (myrank+1)%anzproz, 1, nachbaroben, laenge, MPI_CHAR, (myrank-1+anzproz)%anzproz, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	MPI_Barrier(MPI_COMM_WORLD);	
	
	//printf("kom schwarz in %d\n", myrank);	
	//Untergitter durchgehen, an jedem Punkt Metropolis-Update
	//dann weisse Punkte
	for(int d1=0;d1<teillaenge;d1+=1){
		for(int d2=(teillaenge*anzproz+d1+1)%2;d2<laenge;d2+=2){
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
	
	//printf("sweep weiss in %d\n", myrank);
	//Dann wieder Update Raender
	//Austausch Raender
	for(int i=0;i<laenge;i+=1){
		sendbufferunten[i]=untergitter[(laenge*(teillaenge-1))+i];
		sendbufferoben[i]=untergitter[i];
	}
	MPI_Sendrecv(sendbufferoben, laenge, MPI_CHAR, (myrank-1+anzproz)%anzproz, 0, nachbarunten, laenge, MPI_CHAR, (myrank+1)%anzproz, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	MPI_Sendrecv(sendbufferunten, laenge, MPI_CHAR, (myrank+1)%anzproz, 1, nachbaroben, laenge, MPI_CHAR, (myrank-1+anzproz)%anzproz, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	
	//printf("kom weiss in %d\n", myrank);
	MPI_Allreduce(ergebnisselokal, ergebnisse, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//H, akz, mag berechnen, indem Ergebnis aus jedem Teilgitter aufaddiert wird
	ergebnisse[0]+=H;//Bis jetzt wurde nur Veraenderung des Hamiltonians gemessen->jetzt auch aktualisierung auf Basis des vorherigen Wert
	if(myrank==0){
		double magnetisierung=(ergebnisse[2]>0?ergebnisse[2]/(double)laenge/(double)laenge:-ergebnisse[2]/(double)laenge/(double)laenge);//Betrag bilden
		double magquadrat=magnetisierung*magnetisierung;
		double magvier=magquadrat*magquadrat;
		fprintf(ausgabedatei, "%e\t%e\t%e\t%e\t%e\n",ergebnisse[1]/(double)laenge/(double)laenge, magnetisierung, magquadrat, magvier, ergebnisse[0]);//akz, mag, magquad, magvier, H
	}
	return ergebnisse[0];//Hamiltonian, fuer naechste Rechnung benoetit, zurueckgeben
			
}

void messenmpi(int messungen, int laenge, double T, double j, char *gitter, FILE *ausgabedatei, gsl_rng **generatoren){
	//Fuehrt messungen sweeps am Gitter durch, schreibt Ergebnisse in ausgabedatei
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
	char nachbaroben[laenge];
	char nachbarunten[laenge];
	char untergitter[laenge*teillaenge];		
	//eigenes Untergitter wird mit aktuellem Zustand des Gitters initialisiert. Nur einmal benoetigt, da danach Aenderungen von vorherigen sweep gespeichert sind.
	for (int i=0; i<teillaenge*laenge; i+=1){
		untergitter[i]=gitter[myrank*teillaenge*laenge+i];
	}
	//Raender der benachbarten Untergitter werden gespeichert
	for (int i=0; i<laenge; i+=1){
		nachbaroben[i]=gitter[((teillaenge*myrank-1+laenge)%laenge)*laenge+i];
		nachbarunten[i]=gitter[((teillaenge*(myrank+1))%laenge)*laenge+i];
	}
	//printf("Gitter verteilt in Prozess %d\n", myrank);
	for (int messung=0; messung<messungen; messung+=1){
		if(myrank==0){
			fprintf(ausgabedatei,"%e\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		}
		
		//printf("vor sweep in  %d\n", myrank);
		H=sweepmpi(laenge,ausgabedatei, generatoren,  anzproz, myrank, teillaenge, untergitter, nachbarunten, nachbaroben,  wahrscheinlichkeiten,  j,  H);
	}
	//Gitter aktualisieren nach sweep: In jedem Prozess Inhalte der einzelnen Untergitter in Gitter aneinanderreihen
	MPI_Allgather(untergitter, teillaenge, MPI_CHAR, gitter, teillaenge, MPI_CHAR, MPI_COMM_WORLD);
		
}

void thermalisierenmpi(int messungen, int laenge, double T, double j, char *gitter, FILE *dummymessung, FILE *ausgabedatei, gsl_rng **generatoren){
	//Fuehrt messungen sweeps am Gitter durch, schreibt Ergebnisse in dummymessung, und gibt am Ende Gitter in ausgabedatei
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
	//printf("a %d\n", myrank);
	int teillaenge=laenge/anzproz;
	char nachbaroben[laenge];
	//printf("d %d\n", myrank);
	char nachbarunten[laenge];
	//printf("e %d\n", teillaenge);
	char untergitter[laenge*teillaenge];	
	//printf("b %d\n", myrank);	
	//eigenes Untergitter wird mit aktuellem Zustand des Gitters initialisiert. Nur einmal benoetigt, da danach Aenderungen von vorherigen sweep gespeichert sind.
	for (int i=0; i<teillaenge*laenge; i+=1){
		untergitter[i]=gitter[myrank*teillaenge*laenge+i];
	}
	//printf("c %d\n", myrank);
	//Raender der benachbarten Untergitter werden gespeichert
	for (int i=0; i<laenge; i+=1){
		nachbaroben[i]=gitter[((teillaenge*myrank-1+laenge)%laenge)*laenge+i];
		nachbarunten[i]=gitter[((teillaenge*(myrank+1))%laenge)*laenge+i];
	}
	//printf("Gitter verteilt in Prozess %d\n", myrank);
	for (int messung=0; messung<messungen; messung+=1){
		if(myrank==0){
			fprintf(dummymessung,"%e\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		}		
		//printf("vor funktion in %d\n", myrank);
		H=sweepmpi(laenge,dummymessung, generatoren,  anzproz, myrank, teillaenge, untergitter, nachbarunten, nachbaroben,  wahrscheinlichkeiten,  j,  H);
		//Gitter aktualisieren nach sweep: In jedem Prozess Inhalte der einzelnen Untergitter in Gitter aneinanderreihen
	}
	MPI_Allgather(untergitter, teillaenge, MPI_CHAR, gitter, teillaenge, MPI_CHAR, MPI_COMM_WORLD);
	if(myrank==0){
		ausgabe(gitter, laenge, ausgabedatei);
	}
		
}

void blocks_generieren(int l, int messungen, const int spalte, const int spalten, double *blockarray,  FILE *messdatei){
	//Blockt die Daten aus spalte in Messdatei in Blöcke der Laenge l, Ergebnis in blockarray und ausgabedatei, messdatei muss spalten mit nur doubles enthalten
	double zwischensumme, einwert;//speichern der Summe über die Messwerte und der einzelnen Werte
	double ergebnisarray[spalten];//speichern der ganzen Zeile aus der Messdatei
	rewind(messdatei);
	for (int block=0; block<messungen/l; block+=1){//jedes einzelne Element in Blockarray durchgehen
		zwischensumme=0;
		for (int wert=0; wert<l; wert+=1){//generiert einzelnes Element des blocks
			for (int i=0; i<spalten; i+=1){
				fscanf(messdatei, "%le", &ergebnisarray[i]);
				if (i==spalten-1){fscanf(messdatei, "\n");}
			}
			einwert=ergebnisarray[spalte];//wählt korrekte messung aus
			zwischensumme+=einwert;
		}
		blockarray[block]=zwischensumme/(double)l;
	}

}

double bootstrap_replication(int l, int messungen, double *blockarray, gsl_rng *generator){
	//zieht messungen/l zufaellige Elemente aus blockarray und berechnet den Mittelwert darüber, gibt Mittelwert zurueck
	double summe=0;//Speichert Zwischenergebnis
	int element;////Zufaelliger Integer, der element aus blockarray auswaehlt
	int elementanzahl=(int)messungen/l;//Anzahl der Elemente in blockarray
	for (int durchgang=0; durchgang<elementanzahl; durchgang+=1){
		element=gsl_rng_uniform_int(generator, elementanzahl);//Zufallszahl generieren
		summe+=blockarray[element];//Dazu passenden Messwert auswaehlen
	}
	summe/=(int)(elementanzahl);
	return summe;
}

void bootstrapohnepar(int l, int r, int messungen, double temperatur, double *blockarray, gsl_rng *generator, FILE *ausgabedatei){
//berechnet Mittelwert und Varianz aus r gebootstrappten replikas, schreibt es in ausgabedatei
	double mittelwert=0;//speichern zwischenwerte
	double replica;
	double varianz=0;
	double *bootstraparray;
	if((bootstraparray=(double*)malloc(sizeof(double)*r))==NULL){//speichert ausgewählten Daten, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren des arrays fuer die Replika!\n");
		//return (-1);
	}
	for (int durchgang=0; durchgang<r; durchgang+=1){//Zieht r replicas, speichert sie in bootstraparray und berechnet ihren Mittelwert
		replica=bootstrap_replication(l, messungen, blockarray, generator);
		mittelwert+=replica;
		bootstraparray[durchgang]=replica;//speichern fuer Varianzbildung
	}
	mittelwert/=r;//Standardschaetzer
	for (int durchgang=0; durchgang<r; durchgang+=1){//Berechnet Varianz von ausgewählten werten
		varianz+=(bootstraparray[durchgang]-mittelwert)*(bootstraparray[durchgang]-mittelwert);
	}
	varianz=sqrt(varianz/((double)r-1));//Standardschaetzer
	fprintf(ausgabedatei, "2\t%4d\t%d\t%e\t%e\t%e\n", l,r, mittelwert, varianz, temperatur);//Ausgabe
	free(bootstraparray);
}

extern inline double mittelwertarray(double *array, int messungen){
	//Bestimmt Mittelwert eines arrays, das mit doubles gefüllt ist
	double summe=0;//Speichert Summe über 
	for (int messung=0; messung<messungen; messung+=1){//Mittelwert über Messung bilden
		summe+=array[messung];
	}
	return summe/(double)messungen;
}

extern inline double varianzarray(double *array, int messungen, double mittelwert){
	//Bestimmt die Varianz über alle Elemente in einem array, das mit doubles gefuellt ist, um den Mittelwert
	double summe=0;//Speichert Summe über Messungen
	double einwert=0;//Speichert einen ausgelesenen Wert
	for (int messung=0; messung<messungen; messung+=1){//Mittelwert über Messung bilden
		einwert=array[messung];
		summe+=(einwert-mittelwert)*(einwert-mittelwert);
	}
	return sqrt(summe/((double)messungen-1));
}

extern inline double minarray(double *array, int messungen){
	//Bestimmt Minimum aus array mit gegebener Laenge
	double minwert= array[0];//Fange mit nulltem Wert an
	for (int messung=1; messung<messungen; messung+=1){// gehe alle Elemente durch, fange bei erstem Element an, da nulltes der Startwert ist
		if (array[messung]<minwert){//Vergleich
			minwert=array[messung];//Wenn noetig, Update
		}
	}
	return minwert;
}

extern inline double maxarray(double *array, int messungen){
	//Bestimmt maximum aus array mit gegebener Laenge
	double maxwert= array[0];//Fange mit nulltem Wert an
	for (int messung=1; messung<messungen; messung+=1){// gehe alle Elemente durch, fange bei erstem Element an, da nulltes der Startwert ist
		if (array[messung]>maxwert){//Vergleich
			maxwert=array[messung];//Wenn noetig, Update
		}
	}
	return maxwert;
}


//christiane@christiane-VirtualBox2:17:07:41:~/Bachelorarbeit/bachelorarbeit$ mpicc -g -Wall -o mpifunktionen mpifunktionen.c -lgsl -lm -lgslcblas
//christiane@christiane-VirtualBox2:17:08:52:~/Bachelorarbeit/bachelorarbeit$ mpirun -n 10 ./mpifunktionen
