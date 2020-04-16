//Christiane, start 15.04.20
//Erste Implementierung Ising-Modell für Bachelorarbeit


#include <stdio.h>
#include <gsl/gsl_rng.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_RESET   "\x1b[0m"

void initialisierung(int *gitter, int laenge,/* generatortyp,*/ int seed){
	//initialisiert ein laenge*laenge quadratisches Gitter mit Zufallszahlen -1 und 1
	//initialisiere generator mit seed
	gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(generator, seed);
	unsigned long int zufallsspeicher;
	int zufallsauswertung;	
	int npositiv=0;
	int nnegativ=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			//generiere zufallsspeicher
			zufallsspeicher = gsl_rng_uniform_int(generator, 200);//generiert Zufallszahl zwischen 0 und 200-1
			//printf("%lu ",zufallsspeicher);
			//if zufallsspeicher<max/2: zufallsauswertung =-1
			//else zufallsauswertung =1
			if (zufallsspeicher<100){
				zufallsauswertung=-1;
				nnegativ+=1;
			}
			else{
				zufallsauswertung=1;
				npositiv+=1;}
			gitter[laenge*d1+d2]=zufallsauswertung;//zufallszahl -1 oder 1
			//if (d1+d2<3){gitter[laenge*d1+d2]=-1;}
		}
	}
	gsl_rng_free(generator);
	printf("pos %d, neg %d\n", npositiv, nnegativ);
}

void ausgabe(int *gitter, int laenge){
	//gibt ein laenge *laenge quadratisches Gitter auf die Atandardkonsole aus
	for (int d1=0; d1<laenge; d1++){//geht in erster dimension durch
		printf(" ");
		for (int d2=0; d2<laenge; d2++){//geht in zweiter dimension durch
			if (gitter[laenge*d1+d2]==1){
				printf( ANSI_COLOR_RED "+ " ANSI_COLOR_RESET);
			}
			if (gitter[laenge*d1+d2]==-1){
				printf( ANSI_COLOR_GREEN "- " ANSI_COLOR_RESET);
			}
//else {printf ("E ");}
			//printf("%d ",gitter[laenge*d1+d2]); //einfache ausgabe
		}
		printf("\n");//neue Zeile
	}
}

double hamiltonian(int *gitter, int laenge, double j){
	//berechnet den hamiltonian eines laenge*laenge quadratischen Gitters eines Ising Modells mit periodischen Randbedingungen
	double H=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
		//koennte effizienter sein, eine if-else Schleife für den Schritt nach rechts und eine für den Schritt nach unten zu machen
			if(d2==laenge-1){
				H+=gitter[laenge*d1+d2]*gitter[laenge*d1];//Schritt nach Rechts am Rand
				//printf("rechtsrand ");
			}
			else{
				H+=gitter[laenge*d1+d2]*gitter[laenge*d1+(d2+1)];//Schritt nach Rechts
				//printf("rechts ");
				}
			if(d1==laenge-1){
				H+=gitter[laenge*d1+d2]*gitter[d2];//Schritt nach unten am Rand
				//printf("untenrand ");
				}
			else{
				H+=gitter[laenge*d1+d2]*gitter[laenge*(d1+1)+d2];//Schritt nach unten
				//printf("unten ");
				}
				//printf("%d, %d \n", d1, d2);
		}
	}
	return -j*H;
	
}

double deltah(int *gitter, int d1, int d2, int laenge){
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
	
	return delta;
}

void flipspin(int *gitter, int d1, int d2, int laenge){
	gitter[laenge*d1+d2]*=-1;
}
	

int main(int argc, char **argv){
	int laenge=30;
	double j=1.0;
	//typ=;
	//seed=;
	int gitter [laenge*laenge];//initialisiert gitter
	initialisierung(gitter, laenge, 5);//, typ, seed);
	ausgabe(gitter, laenge);
	double H=hamiltonian(gitter, laenge, j);
	printf("H=%f\n", H);
	printf("dH=%f", deltah(gitter, 1,1,laenge));
	flipspin(gitter, 1, 1, laenge);
	H=hamiltonian(gitter, laenge, j);
	printf("H=%f\n", H);
	return 0;
}

//christiane@christiane-VirtualBox2:/media/christiane/BC20-2E26/Bachelorarbeit$ gcc -std=c99 -Wall -pedantic -I /usr/include/ ising.c -c
//christiane@christiane-VirtualBox2:/media/christiane/BC20-2E26/Bachelorarbeit$ gcc -std=c99 -Wall -pedantic -o ising.exe ising.o -lgsl -lgslcblas
//christiane@christiane-VirtualBox2:/media/christiane/BC20-2E26/Bachelorarbeit$ ./ising.exe

//verworfener code
			//~ if (d2==laenge-1&&d1!=laenge-1){//periodische Bedingung rechts
				//~ H+=gitter[laenge*d1+d2]*gitter[laenge*d1];//Schritt nach Rechts
				//~ H+=gitter[laenge*d1+d2]*gitter[laenge*(d1+1)+d2];//Schritt nach unten
			//~ }
			//~ if (d1==laenge-1&&d2!=laenge-1){//periodische Bedingung unten
				//~ H+=gitter[laenge*d1+d2]*gitter[laenge*d1+(d2+1)];//Schritt nach Rechts
				//~ H+=gitter[laenge*d1+d2]*gitter[d2];//Schritt nach unten
			//~ }
			//~ if (d1==laenge-1&&d2==laenge-1){//periodische Bedingung unten rechts
				//~ H+=gitter[laenge*d1+d2]*gitter[laenge*d1];//Schritt nach Rechts
				//~ H+=gitter[laenge*d1+d2]*gitter[d2];//Schritt nach unten
			//~ }
			//~ else{
				//~ H+=gitter[laenge*d1+d2]*gitter[laenge*d1+(d2+1)];//Schritt nach Rechts
				//~ 
				//~ }

