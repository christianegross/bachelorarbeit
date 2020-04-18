//Christiane, start 15.04.20
//Erste Implementierung Ising-Modell für Bachelorarbeit
//Versuch Metropolis-Algorithmus


#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "math.h"

#define ANSI_COLOR_RED     "\x1b[31m" //methode gefunden auf https://stackoverflow.com/questions/3219393/stdlib-and-colored-output-in-c
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RESET   "\x1b[0m"

void initialisierung(int *gitter, int laenge, int seed){
	//initialisiert ein laenge*laenge quadratisches Gitter mit Zufallszahlen -1 und 1
	//initialisiere generator mit seed
	gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	gsl_rng_set(generator, seed);
	unsigned long int zufallsspeicher;
	int zufallsauswertung;
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


void ausgabebunt(int *gitter, int laenge){
	//gibt ein laenge *laenge quadratisches Gitter auf die Standardkonsole aus
	//farbige Ausgabe für +-1
	for (int d1=0; d1<laenge; d1++){//geht in erster dimension durch
		printf(" ");
		for (int d2=0; d2<laenge; d2++){//geht in zweiter dimension durch
			switch(gitter[laenge*d1+d2]){
				case 1:
					printf( ANSI_COLOR_RED "+ " ANSI_COLOR_RESET);
					break;
				case -1:
					printf( ANSI_COLOR_GREEN "- " ANSI_COLOR_RESET);
					break;
				default:
					printf("%d ", gitter[laenge*d1+d2]);//, gitter[laenge*d1+d2]);
					break;
			}
		}
		printf("\n");//neue Zeile
	}
}

void ausgabe(int *gitter, int laenge, FILE *datei){
	//gibt ein laenge *laenge quadratisches Gitter in datei aus
	for (int d1=0; d1<laenge; d1++){//geht in erster dimension durch
		for (int d2=0; d2<laenge; d2++){//geht in zweiter dimension durch
			fprintf(datei, "%d \t %d \t %d \n",d1, d2, gitter[laenge*d1+d2]);//, gitter[laenge*d1+d2]);
			//printf("%d \t %d \t %d \n",d1, d2, gitter[laenge*d1+d2]);//, gitter[laenge*d1+d2]);
		}
	}
}

void einlesen(int *gitter, int laenge, FILE *datei){
	//liest gitter, das in datei geschrieben wurde, wieder in gitter ein
	//int d1[]={0};
	//int d2[]={0};
	rewind(datei);
	int d1, d2, error;
	for (int n=0; n<laenge*laenge; n+=1){
		error= fscanf(datei, "%d \t %d \t %d \n", &d1, &d2, &gitter[n]);
		//printf("%d \t %d \t %d \n", d1, d2,gitter[n]);
		//gitter[n]=//dritter gelesener Eintrag in Zeile
	}
}
		

double hamiltonian(int *gitter, int laenge, double j){
	//berechnet den hamiltonian eines laenge*laenge quadratischen Gitters eines Ising Modells mit periodischen Randbedingungen
	double H=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			H+=gitter[laenge*d1+d2]*gitter[laenge*d1+((d2+1)%(laenge))];
			H+=gitter[laenge*d1+d2]*gitter[laenge*((d1+1)%(laenge))+d2];
		}
	}
	return -j*H;	
}


double deltah(int *gitter, int d1, int d2, int laenge, double j){
	//Nur eine Variable zur Angabe der Position? oben/unten-+laenge, rechts/links+-1, randbestimmung mit modulo
	//berechnet Energieänderung bei Flip des Spins an position d1, d2
	double delta=0;
	//-2*aktueller Zustand: 1-(2*1)=-1, (-1)-(-1*2)=1
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*((d1-1+laenge)%(laenge))+d2];//oben
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*((d1+1)%(laenge))+d2];//unten
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*d1+((d2-1+laenge)%(laenge))];//links
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*d1+((d2+1)%(laenge))];//rechts
	return -j*delta;
}

int tryflip(int *gitter,  int d1, int d2, int laenge, double j, double T, gsl_rng *generator, double delta){
	//versucht, den spin an position d1, d2 umzukehren nach Metropolis-Algorithmus
	//if deltah<0: accept, return 1
	if (delta<0){
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

void flipspin(int *gitter, int d1, int d2, int laenge){
	gitter[laenge*d1+d2]*=-1;
}
	
double sweep(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *fp){
	double H=hamiltonian;
	double delta;
	int changes=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			delta=deltah(gitter, d1, d2, laenge, j);
			if (tryflip(gitter, d1, d2, laenge, j, T, generator, delta)==1){
				flipspin(gitter, d1, d2, laenge);
				H+=delta;
				//printf("H=%f\n", H);
				changes+=1;
			}
		}
	}
	//ausgabe(gitter, laenge);
	//printf("changes=%d of %d possibilities\n", changes, laenge*laenge);
	fprintf(fp, "%d\n", changes);//benoetigte messungen
	return H;
}

//funktion "thermalisieren" benötigt?
//Für große Gitter: thermalisieren dauert lange, Ergebnis in Datei speichern, nachher einlesen

int main(int argc, char **argv){
	int laenge=300;
	double j=1.0;
	int seed=5;
	//double gamma_T=10.0; //benutzt für SA
	//int numberofsweeps, temperaturechanges;
	double T=10;
	//int messungen=1000;
	
	int gitter [laenge*laenge];//initialisiert gitter
	initialisierung(gitter, laenge, seed);
	//ausgabe(gitter, laenge);
	double H=hamiltonian(gitter, laenge, j); //Anfangsenergie
	printf("H=%f\n", H);
	
	
	gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	gsl_rng_set(generator, seed);
	double Hneu=H;
	double Halt=H+laenge*j+1;
	//int sumdifferences;
	//thermalisieren
	FILE *therm=fopen("thermalisierung.txt", "w");
	FILE *gittertherm=fopen("gitterthermalisiert.txt", "w+");
	int N0=0;//zählt, wie viele sweeps zum Thermalisieren benoetigt werden
	while (Halt-Hneu>0){//Abbruchkriterium
		Halt=Hneu;//Zustand der vorherigen Iteration speichern zum Vergleich
		Hneu=sweep(gitter, laenge, j, T, generator, Halt, therm);//neuen Zustand durch sweep vom alten Zustand
		fprintf(therm,"%d %f\n", N0, Halt-Hneu);
		N0+=1;
	}
	ausgabe(gitter, laenge, gittertherm);
	//fclose(gittertherm);
	//gittertherm=fopen("gitterthermalisiert.txt", "r");
	einlesen(gitter, laenge, gittertherm);
	//printf("%d \n", N0);
	fclose(gittertherm);
	fclose(therm);
	//messen
	//~ FILE *messen=fopen("messen.txt", "w");//speichert messungen
	//~ for (int messung=0; messung<messungen; messung+=1){
		//~ fprintf(messen,"%d\t", messung);
		//~ H=sweep(gitter, laenge, j, T, generator, H, messen);
	//~ }	
	//~ fclose(messen);
//Metropolis:
//N_0 sweeps zum thermalisieren
//Wann thermalisiert? Schwelle für Änderung Hamiltonian?
//Plot: Änderung Hamiltonian wird nicht geringer für T=10
//Fuer T kleiner, kurz große Änderungen, aber danach imemr noch "grosse" Schwankungen
//Fuer groessere laengen Effekt sichtbarer, Abbruchkriterium: Diff=0
//N sweeps zum Messen: Was als Rückgabewert von sweep? sweep überhaupt als Funktion oder in main integrieren?
//Als erstes messen: Akzeptanzrate
	
	gsl_rng_free(generator);

	return 0;
}

//Ideen für Ausgabe über Gnuplot: 
//Heatmap machen, brauche x, y Koordinate, Wert
//x=d1, y=d2, Wert=+-1
//schreibe für initialisieren/sweep in Datei, entweder viele eigene oder durch zwei Leerzeilen in Blöcke getrennt
//Code Heatmap wie in Versuch 362
//Ausgabe der Gnuplot commands durch printf in Datei?
//In Datei schreiben: Siehe Computerphysik
	//~ double T=gamma_T; //simulated annealing
	//~ for (temperaturechanges=0; temperaturechanges<20; temperaturechanges+=1){//Temperatur wird zehnmal aktualisiert
		//~ for (numberofsweeps=0;numberofsweeps<10;numberofsweeps+=1){//zehn sweeps pro Temperatur
			//~ H=sweep(gitter, laenge, j, T, generator, H);
		//~ }
		//~ T=gamma_T/log(temperaturechanges+1.5);//+1.5, da log(1)=0, implicit conversion in double?
		//~ printf("H=%f\n", H);
		//~ ausgabe(gitter, laenge);
	//~ }
//läuft nur mit explizitem Suchen nach gsl und anhängen von zusätzlichen bibliotheken
//christiane@christiane-VirtualBox2:/media/christiane/BC20-2E26/Bachelorarbeit$ gcc -std=c99 -Wall -pedantic -I /usr/include/ ising.c -c
//christiane@christiane-VirtualBox2:/media/christiane/BC20-2E26/Bachelorarbeit$ gcc -std=c99 -Wall -pedantic -o ising.exe ising.o -lgsl -lgslcblas -lm
//christiane@christiane-VirtualBox2:/media/christiane/BC20-2E26/Bachelorarbeit$ ./ising.exe
