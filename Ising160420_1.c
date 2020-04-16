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


void ausgabe(int *gitter, int laenge){
	//gibt ein laenge *laenge quadratisches Gitter auf die Standardkonsole aus
	//Gitter darf nur +-1 enthalten
	for (int d1=0; d1<laenge; d1++){//geht in erster dimension durch
		printf(" ");
		for (int d2=0; d2<laenge; d2++){//geht in zweiter dimension durch
			if (gitter[laenge*d1+d2]==1){
				printf( ANSI_COLOR_RED "+ " ANSI_COLOR_RESET);
			}
			if (gitter[laenge*d1+d2]==-1){
				printf( ANSI_COLOR_GREEN "- " ANSI_COLOR_RESET);
			}
		}
		printf("\n");//neue Zeile
	}
}



double hamiltonian(int *gitter, int laenge, double j){
	//berechnet den hamiltonian eines laenge*laenge quadratischen Gitters eines Ising Modells mit periodischen Randbedingungen
	double H=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			if(d2==laenge-1){
				H+=gitter[laenge*d1+d2]*gitter[laenge*d1];//Schritt nach Rechts am Rand
			}
			else{
				H+=gitter[laenge*d1+d2]*gitter[laenge*d1+(d2+1)];//Schritt nach Rechts
				}
			if(d1==laenge-1){
				H+=gitter[laenge*d1+d2]*gitter[d2];//Schritt nach unten am Rand
				}
			else{
				H+=gitter[laenge*d1+d2]*gitter[laenge*(d1+1)+d2];//Schritt nach unten
				}
		}
	}
	return -j*H;	
}

double hamiltonianmod(int *gitter, int laenge, double j){
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
	if (d1==0){
		delta-=2*gitter[laenge*d1+d2]*gitter[laenge*(laenge-1)+d2];//oben Rand
	}
	else{
		delta-=2*gitter[laenge*d1+d2]*gitter[laenge*(d1-1)+d2];//oben
	}
	
	if (d1==laenge-1){
		delta-=2*gitter[laenge*d1+d2]*gitter[d2];//unten Rand
	}
	else{
		delta-=2*gitter[laenge*d1+d2]*gitter[laenge*(d1+1)+d2];//unten
	}
	
	if (d2==0){
		delta-=2*gitter[laenge*d1+d2]*gitter[laenge*d1+(laenge-1)];//links Rand
	}
	else{
		delta-=2*gitter[laenge*d1+d2]*gitter[laenge*d1+(d2-1)];//links
	}
	
	if(d2==laenge-1){
		delta-=2*gitter[laenge*d1+d2]*gitter[laenge*d1];//Rechts Rand
	}
	else{
		delta-=2*gitter[laenge*d1+d2]*gitter[laenge*d1+(d2+1)];//rechts
	}
	
	return -j*delta;
}

int tryflip(int *gitter,  int d1, int d2, int laenge, double j, double T, gsl_rng *generator, double delta){
	//versucht, den spin an position d1, d2 umzukehren nach Metropolis-Algorithmus
	//berechne exp(-DeltaH/T) Rechenleistung sparen?: deltah vorher berechnen, als argument angeben, danach zum updaten des hamiltonians verwenden
	//double delta=deltah(gitter, d1, d2, laenge, j);
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
	
double sweep(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian){
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
	printf("changes=%d\n", changes);
	return H;
}

int main(int argc, char **argv){
	int laenge=30;
	double j=1.0;
	int seed=5;
	double gamma_T=10.0;
	int numberofsweeps, temperaturechanges;
	
	int gitter [laenge*laenge];//initialisiert gitter
	initialisierung(gitter, laenge, seed);
	ausgabe(gitter, laenge);
	double H=hamiltonian(gitter, laenge, j);
	double Hmod=hamiltonianmod(gitter, laenge, j);
	printf("H=%f\n", H);
	printf("Hmod=%f\n", Hmod);
	
	
	//double delta=deltah(gitter, 4,0,laenge, j);
	//printf("dH=%f\n", delta);
	//flipspin(gitter, 1, 1, laenge);
	//H=hamiltonian(gitter, laenge, j);
	//ausgabe(gitter, laenge);
	//printf("H=%f\n", H);
	
	
	gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	gsl_rng_set(generator, seed);
	double T=gamma_T;
	for (temperaturechanges=0; temperaturechanges<10; temperaturechanges+=1){
		for (numberofsweeps=0;numberofsweeps<10;numberofsweeps+=1){
			H=sweep(gitter, laenge, j, T, generator, H);
		}
		T=gamma_T/log(temperaturechanges+1.5);//+1.5, da log(1)=0, implicit conversion in double?
		printf("H=%f\n", H);
		printf("H=%f\n", hamiltonian(gitter, laenge, j));
		printf("Hmod=%f\n", hamiltonianmod(gitter, laenge, j));
		//ausgabe(gitter, laenge);
	}
	//H=sweep(gitter, laenge, j, T, generator, H); //ausprobieren von sweep
	//printf("H=%f\n", H);
	//~ if (tryflip(gitter, 4, 0, laenge, j, T, generator)==1){ //ausprobieren von tryflip
		//~ flipspin(gitter, 4, 0, laenge);
		//~ H+=delta;
		//~ printf("H=%f\n", H);
	//~ }
	gsl_rng_free(generator);
	//gitter[2*laenge+2]=3;//ausprobieren, was ausgabe mit werten !=+-1 macht
	//ausgabe(gitter, laenge);
	return 0;
}

//läuft nur mit explizitem Suchen nach gsl und anhängen von zusätzlichen bibliotheken
//christiane@christiane-VirtualBox2:/media/christiane/BC20-2E26/Bachelorarbeit$ gcc -std=c99 -Wall -pedantic -I /usr/include/ ising.c -c
//christiane@christiane-VirtualBox2:/media/christiane/BC20-2E26/Bachelorarbeit$ gcc -std=c99 -Wall -pedantic -o ising.exe ising.o -lgsl -lgslcblas -lm
//christiane@christiane-VirtualBox2:/media/christiane/BC20-2E26/Bachelorarbeit$ ./ising.exe


//periodische Randbedingungen durch modulo?
//Problem: bei pos=laenge-1 muss null zurückgegeben werden, ansonsten pos+1
//Keine Zuordnung gefunden, die gleichzeitig nach pos=laenge-1 de Schleife verlässt.



//~ double hamiltonianmod(int *gitter, int laenge, double j){
	//~ //berechnet den hamiltonian eines laenge*laenge quadratischen Gitters eines Ising Modells mit periodischen Randbedingungen
	//~ double H=0;
	//~ int pos1, pos2;
	//~ for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		//~ for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			//~ pos1=d1;
			//~ pos2=d2;
			//~ if(d2==laenge-1){
				//~ pos2=-1;//Schritt nach Rechts am Rand
			//~ }
			//~ H+=gitter[laenge*d1+d2]*gitter[laenge*pos1+(pos2+1)];//Schritt nach Rechts
			//~ if(d1==laenge-1){
				//~ pos1=-1;//Schritt nach unten am Rand
				//~ }
			//~ H+=gitter[laenge*d1+d2]*gitter[laenge*(pos1+1)+pos2];//Schritt nach unten
		//~ }
	//~ }
	//~ return -j*H;	
//~ }


//gibt Ausgaben wie diese: 
 //~ - o- o+ - o+ - o
 //~ - o- o+ - o+ - o
 //~ - o- oo- o
 //~ + - o- o+ - o- o

//~ void ausgabe(int *gitter, int laenge){
	//~ //gibt ein laenge *laenge quadratisches Gitter auf die Standardkonsole aus
	//~ //farbige Ausgabe für +-1
	//~ for (int d1=0; d1<laenge; d1++){//geht in erster dimension durch
		//~ printf(" ");
		//~ for (int d2=0; d2<laenge; d2++){//geht in zweiter dimension durch
			//~ switch(gitter[laenge*d1+d2]){
				//~ case 1:
					//~ printf( ANSI_COLOR_RED "+ " ANSI_COLOR_RESET);
				//~ case -1:
					//~ printf( ANSI_COLOR_GREEN "- " ANSI_COLOR_RESET);
				//~ default:
					//~ printf("o");//, gitter[laenge*d1+d2]);
			//~ }
		//~ }
		//~ printf("\n");//neue Zeile
	//~ }
//~ }

