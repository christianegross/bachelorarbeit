//Christiane, ab 25.05.20
//Funktionen fuer die Bachelorarbeit, die zum Messen des Ising-Modells benötigt werden

#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "messfunktionen1606.h"
#include "sweeps1606.h"


void initialisierung(char *gitter, int laenge, int seed){
	//initialisiert ein laenge*laenge quadratisches Gitter mit Zufallszahlen -1 und 1
	//initialisiere generator mit seed
	gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	gsl_rng_set(generator, seed);
	unsigned long int zufallsspeicher;
	char zufallsauswertung;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			zufallsspeicher = gsl_rng_uniform_int(generator, 200);//generiert Zufallszahl zwischen 0 und 200-1
			if (zufallsspeicher<100){//Zuordnung auf +-1
				zufallsauswertung=-1;
			}
			else{
				zufallsauswertung=1;}
			gitter[laenge*d1+d2]=zufallsauswertung;//zufallszahl -1 oder 1
			}
	}
	gsl_rng_free(generator);
}


void ausgabemitplot(char *gitter, int laenge, FILE *datei, FILE *plotdatei){
	//gibt ein laenge *laenge quadratisches Gitter in datei aus
	//Genau wie ausgabe, nur mit Ausgabe zusaetzlich in ints, damit gnuplot plotten kann
	for (int d1=0; d1<laenge; d1++){//geht in erster dimension durch
		for (int d2=0; d2<laenge; d2++){//geht in zweiter dimension durch
			fprintf(datei, "%3d\t%3d\t%c\n",d1, d2, gitter[laenge*d1+d2]);//Gibt Zeile, Spalte und Wert an
			fprintf(plotdatei, "%3d\t%3d\t%d\n",d1, d2, (int)gitter[laenge*d1+d2]);//Gibt Zeile, Spalte und Wert an
			//printf("%3d\t%3d\t%c\n",d1, d2, gitter[laenge*d1+d2]);//Gibt Zeile, Spalte und Wert an
			//printf("%c\t%c\t%c\n", 1, -1, 100);
			//if ((d1+d2)%2==0){printf("0\n");}
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

int gittersummeohnepar (char *gitter, int laenge){
	//berechnet Summe aller Elemente eines Gitter mit laenge*laenge
	int summe=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			summe+=(int)gitter[laenge*d1+d2];
		}
	}
	return abs(summe);
}	
	
int gittersumme (char *gitter, int laenge){
	//berechnet Summe aller Elemente eines Gitter mit laenge*laenge, parallelisiert dabei nach Moeglichkeit
	int summe=0;
	int zwischensumme =0;
	#pragma omp parallel firstprivate (zwischensumme) shared (summe)
	{
		#pragma omp for
		for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
			for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				zwischensumme+=(int)gitter[laenge*d1+d2];
			}
		}
		#pragma omp critical
		{summe+=zwischensumme;}
	}
	return abs(summe);
}		

double hamiltonian(char *gitter, int laenge, double j){
	//berechnet den hamiltonian eines laenge*laenge quadratischen Gitters eines Ising Modells mit periodischen Randbedingungen
	double H=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			H+=(double)gitter[laenge*d1+d2]*gitter[laenge*d1+((d2+1)%(laenge))];//Bond mit rechtem Nachbar
			H+=(double)gitter[laenge*d1+d2]*gitter[laenge*((d1+1)%(laenge))+d2];//Bond mit unterem Nachbar
			//Von allen Feldern rechts und unten berücksichtig, periodische Randbedingungen->alles berücksichtigt
		}
	}
	return -j*H;	
}

double deltahalt(char *gitter, int d1, int d2, int laenge, double j){
	//berechnet Energieänderung bei Flip des Spins an position d1, d2
	double delta=0;
	//-2*aktueller Zustand: 1-(2*1)=-1, (-1)-(-1*2)=1
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*((d1-1+laenge)%(laenge))+d2];//oben
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*((d1+1)%(laenge))+d2];//unten
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*d1+((d2-1+laenge)%(laenge))];//links
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*d1+((d2+1)%(laenge))];//rechts
	return -j*delta;
}

int tryflipalt(double T, gsl_rng *generator, double delta){
	//versucht, den spin an position d1, d2 umzukehren nach Metropolis-Algorithmus
	//if deltah<0: accept, return 1
	if (delta<=0){
		return 1;
	}
	else{
		double probability=exp(-delta/T);
	//Zufallszahl zwischen null und eins
		double random=gsl_rng_uniform(generator);
		if (random<probability){ 
	//if accepted return 1
			return 1;
		}
		else{
	//if rejected return 0
			return 0;
		}
	}
	return -1;
}


int deltah(char *gitter, int d1, int d2, int laenge){
	//berechnet Energieänderung bei Flip des Spins an position d1, d2
	int delta=0;
	//-2*aktueller Zustand: 1-(2*1)=-1, (-1)-(-1*2)=1
	delta+=2*gitter[laenge*d1+d2]*gitter[laenge*((d1-1+laenge)%(laenge))+d2];//oben
	delta+=2*gitter[laenge*d1+d2]*gitter[laenge*((d1+1)%(laenge))+d2];//unten
	delta+=2*gitter[laenge*d1+d2]*gitter[laenge*d1+((d2-1+laenge)%(laenge))];//links
	delta+=2*gitter[laenge*d1+d2]*gitter[laenge*d1+((d2+1)%(laenge))];//rechts
	return delta;
}

int deltahneu(char *gitter, int d1, int d2, int laenge){
	//berechnet Energieänderung bei Flip des Spins an position d1, d2
	int delta=0;
	//-2*aktueller Zustand: 1-(2*1)=-1, (-1)-(-1*2)=1
	char aktuellerwert=gitter[laenge*d1+d2];
	delta+=2*aktuellerwert*gitter[laenge*((d1-1+laenge)%(laenge))+d2];//oben
	delta+=2*aktuellerwert*gitter[laenge*((d1+1)%(laenge))+d2];//unten
	delta+=2*aktuellerwert*gitter[laenge*d1+((d2-1+laenge)%(laenge))];//links
	delta+=2*aktuellerwert*gitter[laenge*d1+((d2+1)%(laenge))];//rechts
	return delta;
}

extern inline int deltahneu2(char *gitter, int d1, int d2, int laenge){
	//berechnet Energieänderung bei Flip des Spins an position d1, d2
	int delta=0;
	//-2*aktueller Zustand: 1-(2*1)=-1, (-1)-(-1*2)=1
	//char aktuellerwert=gitter[laenge*d1+d2];
	delta+=2*gitter[laenge*d1+d2]*(gitter[laenge*((d1-1+laenge)%(laenge))+d2]//oben
								  +gitter[laenge*((d1+1)%(laenge))+d2]//unten
								  +gitter[laenge*d1+((d2-1+laenge)%(laenge))]//links
								  +gitter[laenge*d1+((d2+1)%(laenge))]);//rechts
	return delta;
}

int deltahlookup(char *gitter, int d1, int d2, int laenge, int *lookupplus, int *lookupminus){
	//Lookuptable wie in Binder, Heermann vorgeschlagen
	//macht Berechnungen auf VM schneller, auf qbig allerdings nicht und verringert speedup
	//berechnet Energieänderung bei Flip des Spins an position d1, d2
	int delta=0;
	//-2*aktueller Zustand: 1-(2*1)=-1, (-1)-(-1*2)=1
	//char aktuellerwert=gitter[laenge*d1+d2];
	delta+=2*gitter[laenge*d1+d2]*(gitter[laenge*lookupminus[d1]+d2]//oben
								  +gitter[laenge*lookupplus[d1]+d2]//unten
								  +gitter[laenge*d1+lookupminus[d2]]//links
								  +gitter[laenge*d1+lookupplus[d2]]);//rechts
	return delta;
}

extern inline int tryflip(gsl_rng *generator, double wahrscheinlichkeit){
	//versucht, den spin an position d1, d2 umzukehren nach Metropolis-Algorithmus
	//if deltah<0: accept, return 1
	//~ if (wahrscheinlichkeit==1){
		//~ return 1;
	//~ }
	//~ else{
	//Zufallszahl zwischen null und eins
		double random=gsl_rng_uniform(generator);
		if (random<wahrscheinlichkeit){ 
	//if accepted return 1
			return 1;
		}
		else{
	//if rejected return 0
			return 0;
		}

	//~ }
	return -1;
}

void flipspin(char *gitter, int d1, int d2, int laenge){
	gitter[laenge*d1+d2]*=-1;
}

double wahrscheinlichkeit(int delta, double *wahrscheinlichkeiten){
	double wahrscheinlichkeitwert;
	switch(delta){
	case -8:
		wahrscheinlichkeitwert=wahrscheinlichkeiten[0];
		break;
	case -4:
		wahrscheinlichkeitwert=wahrscheinlichkeiten[1];
		break;
	case 0:
		wahrscheinlichkeitwert=wahrscheinlichkeiten[2];
		break;
	case 4:
		wahrscheinlichkeitwert=wahrscheinlichkeiten[3];
		break;
	case 8:
		wahrscheinlichkeitwert=wahrscheinlichkeiten[4];
		break;
	default: 
		printf("Fehler bei der Berechnung von Delta!\t%d\n", delta);
		wahrscheinlichkeitwert=-1;
		break;
	}
	return wahrscheinlichkeitwert;
}
	


void thermalisierenmitplot(int laenge, double T, double j, int seed,int N0, char *gitter, FILE *ausgabedatei, FILE *plotdatei, gsl_rng *generator){
	//erzeugt ein thermalisiertes Gitter mit laenge*laenge, T, j, seed in ausgabedatei
	//seed im Moment nicht benoetigt, da Gitter von vorheriger Temperatur benutzt wird
	//Genau wie thermalisieren, nur mit Ausgabe zusaetzlich in ints, damit gnuplot plotten kann
	//generator für thermalisieren innerhalb derFunktion seeden?
	//int gitter[laenge*laenge];
	//initialisierung(gitter, laenge, seed);//Initialisiert Gitter
	double H=hamiltonian(gitter, laenge, seed);//Anfangsenergie
	double Hneu=H;
	double Halt=H+laenge*j+1;
	FILE *dummyfile=fopen("dummy.txt", "w");//speichert messergebnisse waehrend des thermalisierens->Nicht benötigt
	for (int anzahl=0; anzahl<N0; anzahl+=1){//Thermalisierungskriterium: feste anzhl an sweeps, durch Parameter übergeben
		Halt=Hneu;//Zustand der vorherigen Iteration speichern zum Vergleich
		Hneu=sweep(gitter, laenge, j, T, generator, Halt, dummyfile);//neuen Zustand durch sweep vom alten Zustand
	}
	//printf("%f\t%d\n", T, N0);zum darstellen Schritte gegen Temperatur
	fclose(dummyfile);
	ausgabemitplot(gitter, laenge, ausgabedatei, plotdatei);//Gitter muss nicht immer neu thermalisiert werden, sondern kann auch eingelesen werden
}

void thermalisieren(int laenge, double T, double j, int seed,int N0, char *gitter, FILE *ausgabedatei, gsl_rng *generator){
	//erzeugt ein thermalisiertes Gitter mit laenge*laenge, T, j, seed in ausgabedatei
	//seed im Moment nicht benoetigt, da Gitter von vorheriger Temperatur benutzt wird
	//generator für thermalisieren innerhalb derFunktion seeden?
	//int gitter[laenge*laenge];
	//initialisierung(gitter, laenge, seed);//Initialisiert Gitter
	double H=hamiltonian(gitter, laenge, seed);//Anfangsenergie
	double Hneu=H;
	double Halt=H+laenge*j+1;
	FILE *dummyfile=fopen("dummy.txt", "w");//speichert messergebnisse waehrend des thermalisierens->Nicht benötigt
	for (int anzahl=0; anzahl<N0; anzahl+=1){//Thermalisierungskriterium: feste anzhl an sweeps, durch Parameter übergeben
		Halt=Hneu;//Zustand der vorherigen Iteration speichern zum Vergleich
		Hneu=sweep(gitter, laenge, j, T, generator, Halt, dummyfile);//neuen Zustand durch sweep vom alten Zustand
	}
	//printf("%f\t%d\n", T, N0);zum darstellen Schritte gegen Temperatur
	fclose(dummyfile);
	ausgabe(gitter, laenge, ausgabedatei);//Gitter muss nicht immer neu thermalisiert werden, sondern kann auch eingelesen werden
}

void messen(int laenge, double T, double j, int messungen, char* gitter/*, FILE *gitterdatei*/, FILE *messdatei, gsl_rng *generator){
	//Führt  messungen Messungen an Gitter in gitterdatei durch mit T, j, generator, speichert das Ergebnis in messdatei
	//char gitter[laenge*laenge];
	//einlesen(gitter, laenge, gitterdatei);
	//~ int lookupplus[laenge], lookupminus[laenge];
	//~ #pragma omp parallel for
	//~ for (int element=0;element<laenge;element+=1){
		//~ lookupplus[element]=element+1;
		//~ lookupminus[element]=element-1;
	//~ }
	//~ lookupplus[laenge-1]=0;
	//~ lookupminus[0]=laenge-1;

	double H=hamiltonian(gitter, laenge, j);
	for (int messung=0; messung<messungen; messung+=1){
		fprintf(messdatei,"%f\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		H=sweep(gitter, laenge, j, T, generator, H, messdatei/*, lookupplus, lookupminus*/);//Geht Gitter durch und schreibt Messwerte in Datei
	}
}

void messenmehreregeneratoren(int laenge, double T, double j, int messungen, char* gitter/*, FILE *gitterdatei*/, FILE *messdatei, gsl_rng **generatoren){
	//Führt  messungen Messungen an Gitter in gitterdatei durch mit T, j, generator, speichert das Ergebnis in messdatei
	//char gitter[laenge*laenge];
	//einlesen(gitter, laenge, gitterdatei);
	double H=hamiltonian(gitter, laenge, j);
	for (int messung=0; messung<messungen; messung+=1){
		fprintf(messdatei,"%f\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		H=sweepeineschleife(gitter, laenge, j, T, generatoren, H, messdatei);//Geht Gitter durch und schreibt Messwerte in Datei
	}
}

void messenvergleichen(int laenge, double T, double j, int messungen, FILE *gitterdatei, FILE *messdatei, FILE *vergleichsdatei, gsl_rng *generator){
	//Führt  messungen Messungen an Gitter in gitterdatei durch mit T, j, generator, speichert das Ergebnis in messdatei
	char gitter1[laenge*laenge];
	einlesen(gitter1, laenge, gitterdatei);
	double H1=hamiltonian(gitter1, laenge, j);
	char gitter2[laenge*laenge];
	einlesen(gitter2, laenge, gitterdatei);
	double H2=hamiltonian(gitter2, laenge, j);
	for (int messung=0; messung<messungen; messung+=1){
		fprintf(messdatei,"%f\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		H1=sweepzweipar(gitter1, laenge, j, T, generator, H1, messdatei);//Geht Gitter durch und schreibt Messwerte in Datei
		H2=sweep(gitter2, laenge, j, T, generator, H2, messdatei);//Geht Gitter durch und schreibt Messwerte in Datei
		fprintf(vergleichsdatei, "%f\t%f\t%f\t%f\n", (double)messung, H1, H2, H1-H2);
	}
	
}
