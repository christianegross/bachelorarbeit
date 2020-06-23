//Christiane, ab 16.06.20
//verschiedene Versionen von sweep-Funktionen

#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen
#include "messfunktionen.h"
#include "sweeps.h"

double sweepaltohnepar(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen){
	//geht das ganze Gitter durch und versucht, jeden Spin umzudrehen. Zählt die Veränderungen, misst Akzeptanzrate und Magneitiserung und gibt aktuellen Hamiltonian zurück
	double H=hamiltonian;
	double delta;
	int changes=0;//Zählt, wie oft geflippt wurde
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			delta=deltahalt(gitter, d1, d2, laenge, j);
			if (tryflipalt(T, generator, delta)==1){//Wenn Spin geflippt wurde
				flipspin(gitter, d1, d2, laenge);//in Gitter speichern
				H+=delta;//H aktualisieren
				changes+=1;
			}
		}
	}
	double akzeptanzrate=(double)changes/(double)laenge/(double)laenge;
	double magnetisierung=(double)gittersumme(gitter, laenge)/(double)laenge/(double)laenge;
	fprintf(dateimessungen, "%f\t%f\t%f\t%f\t%f\n",akzeptanzrate, magnetisierung, magnetisierung*magnetisierung, magnetisierung*magnetisierung*magnetisierung*magnetisierung, H  );//benoetigte messungen: Anzahl Veränderungen+Akzeptanzrate=Veränderungen/Möglichkeiten+Magnetisierung
	return H;
}

double sweepalt(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen){
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
			if (((d1+d2)%2==0)&&(tryflipalt(T, generator, delta)==1)){//Wenn schwarzer Punkt und Spin geflippt wurde
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
			if (((d1+d2)%2==1)&&(tryflipalt(T, generator, delta)==1)){//Wenn weisser Punkt und Spin geflippt wurde
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

double sweepzweipar(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen){
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
				if (((d1+d2)%2==0)&&(tryflip(generator, wahrscheinlichkeit(delta, wahrscheinlichkeiten))==1)){//Wenn schwarzer Punkt und Spin geflippt wurde
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
				if((d1+d2)%2==1){
				delta=deltah(gitter, d1, d2, laenge);
				if (j*(double)delta!=deltahalt(gitter, d1, d2, laenge, j)){
					printf("weiß    Fehler bei delta\n");
				}
				if (((d1+d2)%2==1)&&(tryflip(generator, wahrscheinlichkeit(delta, wahrscheinlichkeiten))==1)){//Wenn weisser Punkt und Spin geflippt wurde
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

double sweep(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen){
	//geht erst alle schwarzen und dann alle weissen Punkte des Gitters durch, macht ein Metropolis-Update an jedem Punkt, schreibt Akzeptanzrate und MAgnetisierung in dateimessungen
	//arbeitet parallel in schleifen ueber die einzelnen Farben
	double H=hamiltonian;//misst Gesamtveraenderung
	double veraenderungH=0;//misst Veraenderung in einem parallen Thread
	int delta=0;
	int d1=0, d2=0;
	double wahrscheinlichkeiten[5]={1,1,1,exp(-4*j/T), exp(-8*j/T)};
	if (j<0){
		wahrscheinlichkeiten[1]=wahrscheinlichkeiten[3];
		wahrscheinlichkeiten[0]=wahrscheinlichkeiten[4];
		wahrscheinlichkeiten[3]=1;
		wahrscheinlichkeiten[4]=1;
		}
	int changes =0;//misst Gesamtzahl der spinflips
	int changesklein=0;//misst Spinflips in parallelen Thread
	//int chunk=2;
	//schwarz: d1+d2 gerade
	//int chunksize=(int)ceil((double)laenge/2.0/(double)omp_get_num_threads());
	#pragma omp parallel firstprivate (delta, veraenderungH, changesklein, d1, d2, generator)// shared(H, changes)
	{
		#pragma omp for nowait schedule (static) //Versuche overhead zu reduzieren
		for (d1=0; d1<laenge;d1+=1){
			for (d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				if((d1+d2)%2==0){
				delta=deltahneu2(gitter, d1, d2, laenge);
				//~ if (j*(double)delta!=deltahalt(gitter, d1, d2, laenge, j)){
					//~ printf("schwarz Fehler bei delta\n");
				//~ }
				if (/*((d1+d2)%2==0)&&*/(tryflip(generator, wahrscheinlichkeit(delta, wahrscheinlichkeiten))==1)){//Wenn schwarzer Punkt und Spin geflippt wurde
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
				delta=deltahneu2(gitter, d1, d2, laenge);
				//~ if (j*(double)delta!=deltahalt(gitter, d1, d2, laenge, j)){
					//~ printf("weiß    Fehler bei delta %d %d\n", d1, d2);
				//~ }
				if (/*((d1+d2)%2==1)&&*/(tryflip(generator, wahrscheinlichkeit(delta, wahrscheinlichkeiten))==1)){//Wenn weisser Punkt und Spin geflippt wurde
					//flipspin(gitter, d1, d2, laenge);//in Gitter speichern
					gitter[laenge*d1+d2]*=-1;
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

double sweeplookup(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen,int *lookupplus, int *lookupminus){
	//geht erst alle schwarzen und dann alle weissen Punkte des Gitters durch, macht ein Metropolis-Update an jedem Punkt, schreibt Akzeptanzrate und MAgnetisierung in dateimessungen
	//arbeitet parallel in schleifen ueber die einzelnen Farben
	double H=hamiltonian;//misst Gesamtveraenderung
	double veraenderungH=0;//misst Veraenderung in einem parallen Thread
	int delta=0;
	int d1=0, d2=0;
	double wahrscheinlichkeiten[5]={1,1,1,exp(-4*j/T), exp(-8*j/T)};
	int changes =0;//misst Gesamtzahl der spinflips
	int changesklein=0;//misst Spinflips in parallelen Thread
	//int chunk=2;
	//schwarz: d1+d2 gerade
	//int chunksize=(int)ceil((double)laenge/2.0/(double)omp_get_num_threads());
	#pragma omp parallel firstprivate (delta, veraenderungH, changesklein, d1, d2)// shared(H, changes)
	{
		#pragma omp for nowait schedule (static) //Versuche overhead zu reduzieren
		for (d1=0; d1<laenge;d1+=1){
			for (d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				if((d1+d2)%2==0){
				delta=deltahlookup(gitter, d1, d2, laenge, lookupplus, lookupminus);
				//~ if (j*(double)delta!=deltahalt(gitter, d1, d2, laenge, j)){
					//~ printf("schwarz Fehler bei delta\n");
				//~ }
				if (/*((d1+d2)%2==0)&&*/(tryflip(generator, wahrscheinlichkeit(delta, wahrscheinlichkeiten))==1)){//Wenn schwarzer Punkt und Spin geflippt wurde
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
				delta=deltahlookup(gitter, d1, d2, laenge, lookupplus, lookupminus);
				//~ if (j*(double)delta!=deltahalt(gitter, d1, d2, laenge, j)){
					//~ printf("weiß    Fehler bei delta %d %d\n", d1, d2);
				//~ }
				if (/*((d1+d2)%2==1)&&*/(tryflip(generator, wahrscheinlichkeit(delta, wahrscheinlichkeiten))==1)){//Wenn weisser Punkt und Spin geflippt wurde
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
		

double sweepmehreregeneratoren(char *gitter, int laenge, double j, double T, gsl_rng **generatoren, double hamiltonian, double *wahrscheinlichkeiten, FILE *dateimessungen){
	//geht erst alle schwarzen und dann alle weissen Punkte des Gitters durch, macht ein Metropolis-Update an jedem Punkt, schreibt Akzeptanzrate und MAgnetisierung in dateimessungen
	//arbeitet parallel in schleifen ueber die einzelnen Farben
	//Fuer jeden Thread einen einzelnen Generator
	//initialisierung der inneren Schleifen mit modulo, um leere Durchläufe zu verhindern
	double H=hamiltonian;//misst Gesamtveraenderung
	double veraenderungH=0;//misst Veraenderung in einem parallen Thread
	int delta=0;
	int d1=0, d2=0;
	int changes =0;//misst Gesamtzahl der spinflips
	int changesklein=0;//misst Spinflips in parallelen Thread
	//int chunk=2;
	//schwarz: d1+d2 gerade
	//int chunksize=(int)ceil((double)laenge/2.0/(double)omp_get_num_threads());
	#pragma omp parallel \
	firstprivate (delta, veraenderungH, changesklein, d1, d2, wahrscheinlichkeiten)//wahrscheinlichkeiten[0], wahrscheinlichkeiten[1], wahrscheinlichkeiten[2], wahrscheinlichkeiten[3], wahrscheinlichkeiten[4])// shared(H, changes)
	{
		int threadnummer=omp_get_thread_num();
		#pragma omp for nowait schedule (static) //Versuche overhead zu reduzieren
		for (d1=0; d1<laenge;d1+=1){
			for (d2=(d1%2); d2<laenge; d2+=2){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				delta=deltahneu2(gitter, d1, d2, laenge);
				if ((tryflip(generatoren[threadnummer], wahrscheinlichkeiten[(delta/4)+2])==1)){//Wenn schwarzer Punkt und Spin geflippt wurde
					gitter[laenge*d1+d2]*=-1;
					veraenderungH+=j*delta;//Zwischenvariable, damit es keine Konflikte beim updaten gibt
					changesklein+=1;
				}
			}
		}
		#pragma omp barrier//damit mit nowait overhead reduziert werden kann
		#pragma omp for nowait schedule (static)
		for (d1=0; d1<laenge;d1+=1){
			for (d2=(d1+1)%2; d2<laenge; d2+=2){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				delta=deltahneu2(gitter, d1, d2, laenge);
				if ((tryflip(generatoren[threadnummer], wahrscheinlichkeiten[(delta/4)+2])==1)){//Wenn weisser Punkt und Spin geflippt wurde
					gitter[laenge*d1+d2]*=-1;
					veraenderungH+=j*delta;
					changesklein+=1;
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
	double magquadrat=magnetisierung*magnetisierung;
	double magvier=magquadrat*magquadrat;
	fprintf(dateimessungen, "%f\t%f\t%f\t%f\t%f\n",akzeptanzrate, magnetisierung, magquadrat, magvier, H  );//benoetigte messungen: Anzahl Veränderungen+Akzeptanzrate=Veränderungen/Möglichkeiten+Magnetisierung
	return H;
}


double sweepeineschleife(char *gitter, int laenge, double j, double T, gsl_rng **generatoren, double hamiltonian, FILE *dateimessungen){
	//geht erst alle schwarzen und dann alle weissen Punkte des Gitters durch, macht ein Metropolis-Update an jedem Punkt, schreibt Akzeptanzrate und MAgnetisierung in dateimessungen
	//arbeitet parallel in schleifen ueber die einzelnen Farben
	//Fuer jeden Thread einen einzelnen Generator
	//Eine Schleife pro Farbe, Zuordnung von d1, d2 aus Schleifenvariable
	double H=hamiltonian;//misst Gesamtveraenderung
	double veraenderungH=0;//misst Veraenderung in einem parallen Thread
	int delta=0;
	int d1=0, d2=0;
	double wahrscheinlichkeiten[5]={1,1,1,exp(-4*j/T), exp(-8*j/T)};
	if (j<0){
		wahrscheinlichkeiten[1]=wahrscheinlichkeiten[3];
		wahrscheinlichkeiten[0]=wahrscheinlichkeiten[4];
		wahrscheinlichkeiten[3]=1;
		wahrscheinlichkeiten[4]=1;
		}
	int changes =0;//misst Gesamtzahl der spinflips
	int changesklein=0;//misst Spinflips in parallelen Thread
	int lhalbe=laenge/2;
	//int chunk=2;
	//schwarz: d1+d2 gerade
	//int chunksize=(int)ceil((double)laenge/2.0/(double)omp_get_num_threads());
	#pragma omp parallel firstprivate (delta, veraenderungH, changesklein, d1, d2)// shared(H, changes)
	{
		int threadnummer=omp_get_thread_num();
		#pragma omp for nowait schedule (static) //Versuche overhead zu reduzieren
		for (int pos=0; pos<laenge*lhalbe;pos+=1){
			d1=(pos-(pos%lhalbe))/lhalbe;
			d2=2*(pos-(d1*lhalbe))+(d1%2);
			delta=deltahneu2(gitter, d1, d2, laenge);
			if ((tryflip(generatoren[threadnummer], wahrscheinlichkeiten[delta/4+2])==1)){//Wenn schwarzer Punkt und Spin geflippt wurde
				gitter[laenge*d1+d2]*=-1;
				veraenderungH+=j*delta;//Zwischenvariable, damit es keine Konflikte beim updaten gibt
				changesklein+=1;
			}
		}
		#pragma omp barrier//damit mit nowait overhead reduziert werden kann
		#pragma omp for nowait schedule (static)
		for (int pos=0; pos<laenge*lhalbe;pos+=1){
			d1=(pos-(pos%lhalbe))/lhalbe;
			d2=2*(pos-(d1*lhalbe))+((d1+1)%2);
			delta=deltahneu2(gitter, d1, d2, laenge);
			if ((tryflip(generatoren[threadnummer], wahrscheinlichkeiten[delta/4+2])==1)){//Wenn weisser Punkt und Spin geflippt wurde
				gitter[laenge*d1+d2]*=-1;
				veraenderungH+=j*delta;
				changesklein+=1;
			}
		}
		#pragma omp critical (hamiltonian)//zwei kritische Regionen, damit weniger Zeit geblockt wird
		{H+=veraenderungH;}
		#pragma omp critical (changes)
		{changes+=changesklein;}
		#pragma omp barrier
	}
	double akzeptanzrate=(double)changes/(double)laenge/(double)laenge;
	double magnetisierung=(double)gittersummeohnepar(gitter, laenge)/(double)laenge/(double)laenge;
	fprintf(dateimessungen, "%f\t%f\n",akzeptanzrate, magnetisierung );//benoetigte messungen: Anzahl Veränderungen+Akzeptanzrate=Veränderungen/Möglichkeiten+Magnetisierung
	return H;
}


