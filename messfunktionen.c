//Christiane, ab 25.05.20
//Funktionen fuer die Bachelorarbeit, die zum Messen des Ising-Modells benötigt werden

#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "messfunktionen.h"


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


void ausgabe(int *gitter, int laenge, FILE *datei){
	//gibt ein laenge *laenge quadratisches Gitter in datei aus
	for (int d1=0; d1<laenge; d1++){//geht in erster dimension durch
		for (int d2=0; d2<laenge; d2++){//geht in zweiter dimension durch
			fprintf(datei, "%3d \t %3d \t %d \n",d1, d2, gitter[laenge*d1+d2]);//Gibt Zeile, Spalte und Wert an
		}
	}
}

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

int gittersummeohnepar (int *gitter, int laenge){
	//berechnet Summe aller Elemente eines Gitter mit laenge*laenge
	int summe=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			summe+=gitter[laenge*d1+d2];
		}
	}
	return abs(summe);
}	
	
int gittersumme (int *gitter, int laenge){
	//berechnet Summe aller Elemente eines Gitter mit laenge*laenge, parallelisiert dabei nach Moeglichkeit
	int summe=0;
	int zwischensumme =0;
	#pragma omp parallel firstprivate (zwischensumme) shared (summe)
	{
		#pragma omp for
		for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
			for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				zwischensumme+=gitter[laenge*d1+d2];
			}
		}
		#pragma omp critical
		{summe+=zwischensumme;}
	}
	return abs(summe);
}		

double hamiltonian(int *gitter, int laenge, double j){
	//berechnet den hamiltonian eines laenge*laenge quadratischen Gitters eines Ising Modells mit periodischen Randbedingungen
	double H=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			H+=gitter[laenge*d1+d2]*gitter[laenge*d1+((d2+1)%(laenge))];//Bond mit rechtem Nachbar
			H+=gitter[laenge*d1+d2]*gitter[laenge*((d1+1)%(laenge))+d2];//Bond mit unterem Nachbar
			//Von allen Feldern rechts und unten berücksichtig, periodische Randbedingungen->alles berücksichtigt
		}
	}
	return -j*H;	
}

double deltahalt(int *gitter, int d1, int d2, int laenge, double j){
	//berechnet Energieänderung bei Flip des Spins an position d1, d2
	double delta=0;
	//-2*aktueller Zustand: 1-(2*1)=-1, (-1)-(-1*2)=1
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*((d1-1+laenge)%(laenge))+d2];//oben
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*((d1+1)%(laenge))+d2];//unten
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*d1+((d2-1+laenge)%(laenge))];//links
	delta-=2*gitter[laenge*d1+d2]*gitter[laenge*d1+((d2+1)%(laenge))];//rechts
	return -j*delta;
}

int tryflipalt(int *gitter,  int d1, int d2, int laenge, double j, double T, gsl_rng *generator, double delta){
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


int deltah(int *gitter, int d1, int d2, int laenge){
	//berechnet Energieänderung bei Flip des Spins an position d1, d2
	double delta=0;
	//-2*aktueller Zustand: 1-(2*1)=-1, (-1)-(-1*2)=1
	delta+=2*gitter[laenge*d1+d2]*gitter[laenge*((d1-1+laenge)%(laenge))+d2];//oben
	delta+=2*gitter[laenge*d1+d2]*gitter[laenge*((d1+1)%(laenge))+d2];//unten
	delta+=2*gitter[laenge*d1+d2]*gitter[laenge*d1+((d2-1+laenge)%(laenge))];//links
	delta+=2*gitter[laenge*d1+d2]*gitter[laenge*d1+((d2+1)%(laenge))];//rechts
	return delta;
}

int tryflip(int *gitter,  int d1, int d2, int laenge, double j, double T, gsl_rng *generator, double wahrscheinlichkeit){
	//versucht, den spin an position d1, d2 umzukehren nach Metropolis-Algorithmus
	//if deltah<0: accept, return 1
	if (wahrscheinlichkeit==1){
		return 1;
	}
	else{
		//double probability=exp(-delta/T);
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
	}
	return -1;
}

void flipspin(int *gitter, int d1, int d2, int laenge){
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
		printf("Fehler bei der Berechnung von Delta!\n");
		wahrscheinlichkeitwert=-1;
		break;
	}
	return wahrscheinlichkeitwert;
}
	
double sweepaltohnepar(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen){
	//geht das ganze Gitter durch und versucht, jeden Spin umzudrehen. Zählt die Veränderungen, misst Akzeptanzrate und Magneitiserung und gibt aktuellen Hamiltonian zurück
	double H=hamiltonian;
	double delta;
	int changes=0;//Zählt, wie oft geflippt wurde
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			delta=deltahalt(gitter, d1, d2, laenge, j);
			if (tryflipalt(gitter, d1, d2, laenge, j, T, generator, delta)==1){//Wenn Spin geflippt wurde
				flipspin(gitter, d1, d2, laenge);//in Gitter speichern
				H+=delta;//H aktualisieren
				changes+=1;
			}
		}
	}
	double akzeptanzrate=(double)changes/(double)laenge/(double)laenge;
	double magnetisierung=(double)gittersumme(gitter, laenge)/(double)laenge/(double)laenge;
	fprintf(dateimessungen, "%f\t%f\t",akzeptanzrate, magnetisierung );//benoetigte messungen: Anzahl Veränderungen+Akzeptanzrate=Veränderungen/Möglichkeiten+Magnetisierung
	return H;
}

double sweepalt(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen){
	//geht erst alle schwarzen und dann alle weissen Punkte des Gitters durch, macht ein Metropolis-Update an jedem Punkt, schreibt Akzeptanzrate und MAgnetisierung in dateimessungen
	double H=hamiltonian;
	double veraenderungH=0;
	double delta=0;
	int changes =0;
	//falls parallel: delta private, changes shared
	//schwarz: d1+d2 gerade
	for (int d1=0; d1<laenge;d1+=1){
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			delta=deltahalt(gitter, d1, d2, laenge, j);
			if (((d1+d2)%2==0)&&(tryflipalt(gitter, d1, d2, laenge, j, T, generator, delta)==1)){//Wenn schwarzer Punkt und Spin geflippt wurde
				flipspin(gitter, d1, d2, laenge);//in Gitter speichern
				veraenderungH+=delta;//H aktualisieren
				changes+=1;
			}
		}
	}
	//printf("\nVeraenderung ohne parallel: %f\t", veraenderungH);
	H+=veraenderungH;
	//weiss: d1+d2 ungerade
	for (int d1=0; d1<laenge;d1+=1){
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			delta=deltahalt(gitter, d1, d2, laenge, j);
			if (((d1+d2)%2==1)&&(tryflipalt(gitter, d1, d2, laenge, j, T, generator, delta)==1)){//Wenn weisser Punkt und Spin geflippt wurde
				flipspin(gitter, d1, d2, laenge);//in Gitter speichern
				H+=delta;//H aktualisieren
				changes+=1;
			}
		}
	}
	double akzeptanzrate=(double)changes/(double)laenge/(double)laenge;
	double magnetisierung=(double)gittersumme(gitter, laenge)/(double)laenge/(double)laenge;
	fprintf(dateimessungen, "%f\t%f\t",akzeptanzrate, magnetisierung );//benoetigte messungen: Anzahl Veränderungen+Akzeptanzrate=Veränderungen/Möglichkeiten+Magnetisierung
	return H;
}

double sweepzweipar(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen){
	//geht erst alle schwarzen und dann alle weissen Punkte des Gitters durch, macht ein Metropolis-Update an jedem Punkt, schreibt Akzeptanzrate und MAgnetisierung in dateimessungen
	//arbeitet parallel in schleifen ueber die einzelnen Farben
	double H=hamiltonian;//misst Gesamtveraenderung
	double veraenderungH=0;//misst Veraenderung in einem parallen Thread
	int delta=0;
	int d1=0, d2=0;
	double wahrscheinlichkeiten[5]={1,1,1,exp(-4*j/T), exp(-8*j/T)};
	int changes =0;//misst Gesamtzahl der spinflips
	int changesklein=0;//misst Spinflips in parallelen Thread
	//schwarz: d1+d2 gerade
	#pragma omp parallel firstprivate (delta, veraenderungH, changesklein, d1, d2) shared (H, changes)
	{
		#pragma omp for
		for (d1=0; d1<laenge;d1+=1){
			for (d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				if((d1+d2)%2==0){
				delta=deltah(gitter, d1, d2, laenge);
				if (j*(double)delta!=deltahalt(gitter, d1, d2, laenge, j)){
					printf("schwarz Fehler bei delta\n");
				}
				if (((d1+d2)%2==0)&&(tryflip(gitter, d1, d2, laenge, j, T, generator, wahrscheinlichkeit(delta, wahrscheinlichkeiten))==1)){//Wenn schwarzer Punkt und Spin geflippt wurde
					//flipspin(gitter, d1, d2, laenge);//in Gitter speichern
					gitter[laenge*d1+d2]*=-1;
					veraenderungH+=j*delta;//Zwischenvariable, damit es keine Konflikte beim updaten gibt
					changesklein+=1;
				}
			}
			}
		}
		#pragma omp critical (schwarzepunkte)//damit das updaten keine konflikte verursacht, Name, damit die critical regionen unabhängig voneinander sind 
		{H+=veraenderungH;
			changes+=changesklein;}
		#pragma omp barrier
	}
	veraenderungH=0;//Zuruecksetzen, damit in naechster paralleler Region nur deren Veraenderungen gezaehlt werden
	changesklein=0;
	//weiss: d1+d2 ungerade, sonst analog zu schwarz
	#pragma omp parallel firstprivate (delta, veraenderungH, changesklein, d1, d2) shared (H, changes)
	{
		#pragma omp for
		for (d1=0; d1<laenge;d1+=1){
			for (d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				if((d1+d2)%2==0){
				delta=deltah(gitter, d1, d2, laenge);
				if (j*(double)delta!=deltahalt(gitter, d1, d2, laenge, j)){
					printf("weiß    Fehler bei delta\n");
				}
				if (((d1+d2)%2==1)&&(tryflip(gitter, d1, d2, laenge, j, T, generator, wahrscheinlichkeit(delta, wahrscheinlichkeiten))==1)){//Wenn weisser Punkt und Spin geflippt wurde
					flipspin(gitter, d1, d2, laenge);//in Gitter speichern
					veraenderungH+=j*delta;
					changesklein+=1;
				}
			}
			}
		}
		#pragma omp critical (weissepunkte)
		{H+=veraenderungH;
		changes+=changesklein;}
		#pragma omp barrier
	}
	double akzeptanzrate=(double)changes/(double)laenge/(double)laenge;
	double magnetisierung=(double)gittersummeohnepar(gitter, laenge)/(double)laenge/(double)laenge;
	fprintf(dateimessungen, "%f\t%f\n",akzeptanzrate, magnetisierung );//benoetigte messungen: Anzahl Veränderungen+Akzeptanzrate=Veränderungen/Möglichkeiten+Magnetisierung
	return H;
}

double sweep(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen){
	//geht erst alle schwarzen und dann alle weissen Punkte des Gitters durch, macht ein Metropolis-Update an jedem Punkt, schreibt Akzeptanzrate und MAgnetisierung in dateimessungen
	//arbeitet parallel in schleifen ueber die einzelnen Farben
	double H=hamiltonian;//misst Gesamtveraenderung
	double veraenderungH=0;//misst Veraenderung in einem parallen Thread
	int delta=0;
	int d1=0, d2=0;
	double wahrscheinlichkeiten[5]={1,1,1,exp(-4*j/T), exp(-8*j/T)};
	int changes =0;//misst Gesamtzahl der spinflips
	int changesklein=0;//misst Spinflips in parallelen Thread
	int chunk=2;
	//schwarz: d1+d2 gerade
	//int chunksize=(int)ceil((double)laenge/2.0/(double)omp_get_num_threads());
	#pragma omp parallel firstprivate (delta, veraenderungH, changesklein, d1, d2)
	{
		#pragma omp for nowait schedule (static) //Versuche overhead zu reduzieren
		for (d1=0; d1<laenge;d1+=1){
			for (d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				if((d1+d2)%2==0){
				delta=deltah(gitter, d1, d2, laenge);
				if (j*(double)delta!=deltahalt(gitter, d1, d2, laenge, j)){
					printf("schwarz Fehler bei delta\n");
				}
				if (((d1+d2)%2==0)&&(tryflip(gitter, d1, d2, laenge, j, T, generator, wahrscheinlichkeit(delta, wahrscheinlichkeiten))==1)){//Wenn schwarzer Punkt und Spin geflippt wurde
					//flipspin(gitter, d1, d2, laenge);//in Gitter speichern
					gitter[laenge*d1+d2]*=-1;
					veraenderungH+=j*delta;//Zwischenvariable, damit es keine Konflikte beim updaten gibt
					changesklein+=1;
				}
			}
			}
		}
		//~ #pragma omp critical (schwarzepunkte)//damit das updaten keine konflikte verursacht, Name, damit die critical regionen unabhängig voneinander sind 
		//~ {H+=veraenderungH;
			//~ changes+=changesklein;}
		//~ #pragma omp barrier
	//~ }
	//~ veraenderungH=0;//Zuruecksetzen, damit in naechster paralleler Region nur deren Veraenderungen gezaehlt werden
	//~ changesklein=0;
	//~ //weiss: d1+d2 ungerade, sonst analog zu schwarz
	//~ #pragma omp parallel firstprivate (delta, veraenderungH, changesklein, d1, d2) shared (H, changes)
	//~ {
		#pragma omp barrier//damit mit nowait overhead reduziert werden kann
		#pragma omp for nowait schedule (static)
		for (d1=0; d1<laenge;d1+=1){
			for (d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				if((d1+d2)%2==1){
				delta=deltah(gitter, d1, d2, laenge);
				if (j*(double)delta!=deltahalt(gitter, d1, d2, laenge, j)){
					printf("weiß    Fehler bei delta %d %d\n", d1, d2);
				}
				if (((d1+d2)%2==1)&&(tryflip(gitter, d1, d2, laenge, j, T, generator, wahrscheinlichkeit(delta, wahrscheinlichkeiten))==1)){//Wenn weisser Punkt und Spin geflippt wurde
					flipspin(gitter, d1, d2, laenge);//in Gitter speichern
					veraenderungH+=j*delta;
					changesklein+=1;
				}
			}
			}
		}
		#pragma omp critical (weissepunkte)
		{H+=veraenderungH;
		changes+=changesklein;}
		#pragma omp barrier
	}
	double akzeptanzrate=(double)changes/(double)laenge/(double)laenge;
	double magnetisierung=(double)gittersummeohnepar(gitter, laenge)/(double)laenge/(double)laenge;
	fprintf(dateimessungen, "%f\t%f\n",akzeptanzrate, magnetisierung );//benoetigte messungen: Anzahl Veränderungen+Akzeptanzrate=Veränderungen/Möglichkeiten+Magnetisierung
	return H;
}
		

void thermalisieren(int laenge, double T, double j, int seed,int N0, int *gitter, FILE *ausgabedatei, gsl_rng *generator){
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

void messen(int laenge, double T, double j, int messungen, FILE *gitterdatei, FILE *messdatei, gsl_rng *generator){
	//Führt  messungen Messungen an Gitter in gitterdatei durch mit T, j, generator, speichert das Ergebnis in messdatei
	int gitter[laenge*laenge];
	einlesen(gitter, laenge, gitterdatei);
	double H=hamiltonian(gitter, laenge, j);
	for (int messung=0; messung<messungen; messung+=1){
		fprintf(messdatei,"%f\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		H=sweep(gitter, laenge, j, T, generator, H, messdatei);//Geht Gitter durch und schreibt Messwerte in Datei
	}
}

void messenvergleichen(int laenge, double T, double j, int messungen, FILE *gitterdatei, FILE *messdatei, FILE *vergleichsdatei, gsl_rng *generator){
	//Führt  messungen Messungen an Gitter in gitterdatei durch mit T, j, generator, speichert das Ergebnis in messdatei
	int gitter1[laenge*laenge];
	einlesen(gitter1, laenge, gitterdatei);
	double H1=hamiltonian(gitter1, laenge, j);
	int gitter2[laenge*laenge];
	einlesen(gitter2, laenge, gitterdatei);
	double H2=hamiltonian(gitter2, laenge, j);
	for (int messung=0; messung<messungen; messung+=1){
		fprintf(messdatei,"%f\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		H1=sweep(gitter1, laenge, j, T, generator, H1, messdatei);//Geht Gitter durch und schreibt Messwerte in Datei
		H2=sweep(gitter2, laenge, j, T, generator, H2, messdatei);//Geht Gitter durch und schreibt Messwerte in Datei
		fprintf(vergleichsdatei, "%f\t%f\t%f\t%f\n", (double)messung, H1, H2, H1-H2);
	}
	
}
