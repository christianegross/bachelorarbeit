//Christiane, start 15.04.20
//Erste Implementierung Ising-Modell für Bachelorarbeit
//Versuch Metropolis-Algorithmus


#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "math.h"


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
			fprintf(datei, "%d \t %d \t %d \n",d1, d2, gitter[laenge*d1+d2]);//Gibt Zeile, Spalte und Wert an
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
	
double sweep(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen){
	double H=hamiltonian;
	double delta;
	int changes=0;//Zählt, wie oft geflippt wurde
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			delta=deltah(gitter, d1, d2, laenge, j);
			if (tryflip(gitter, d1, d2, laenge, j, T, generator, delta)==1){//Wenn Spin geflippt wurde
				flipspin(gitter, d1, d2, laenge);//in Gitter speichern
				H+=delta;//H aktualisieren
				changes+=1;
			}
		}
	}
	double akzeptanzrate=(double)changes/(double)laenge/(double)laenge;
	double magnetisierung=(double)gittersumme(gitter, laenge)/(double)laenge/(double)laenge;
	fprintf(dateimessungen, "%f\t%f\n",akzeptanzrate, magnetisierung );//benoetigte messungen: Anzahl Veränderungen+Akzeptanzrate=Veränderungen/Möglichkeiten+Magnetisierung
	return H;
}

//funktion "thermalisieren" benötigt?
//Für große Gitter: thermalisieren dauert lange, Ergebnis in Datei speichern, nachher einlesen

void thermalisieren(int laenge, double T, double j, int seed, FILE *ausgabedatei, gsl_rng *generator){
	//erzeugt ein thermalisiertes Gitter mit laenge*laenge, T, j, seed in ausgabedatei
	//generator für thermalisieren innerhalb derFunktion seeden?
	int gitter[laenge*laenge];
	initialisierung(gitter, laenge, seed);//Initialisiert Gitter
	double H=hamiltonian(gitter, laenge, seed);//Anfangsenergie
	double Hneu=H;
	double Halt=H+laenge*j+1;
	FILE *dummyfile=fopen("dummy.txt", "w");//speichert messergebnisse waehrend des thermalisierens->Nicht benötigt
	int N0=0;//zählt, wie viele sweeps zum Thermalisieren benoetigt werden
	while (Halt-Hneu>0){//Abbruchkriterium
		Halt=Hneu;//Zustand der vorherigen Iteration speichern zum Vergleich
		Hneu=sweep(gitter, laenge, j, T, generator, Halt, dummyfile);//neuen Zustand durch sweep vom alten Zustand
		N0+=1;
	}
	//printf("%f\t%d\n", T, N0);zum darstellen Schritte gegen Temperatur
	fclose(dummyfile);
	ausgabe(gitter, laenge, ausgabedatei);//ermöglicht Ändern des verwendeten Gitter
}

void messen(int laenge, double T, double j, int messungen, FILE *gitterdatei, FILE *messdatei, gsl_rng *generator){
	//Führt  messungen Messungen an Gitter in gitterdatei durch mit T, j, generator, speichert das Ergebnis in messdatei
	//generator für messen innerhalb derFunktion seeden?
	int gitter[laenge*laenge];
	einlesen(gitter, laenge, gitterdatei);
	double H=hamiltonian(gitter, laenge, j);
	for (int messung=0; messung<messungen; messung+=1){
		fprintf(messdatei,"%f\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		H=sweep(gitter, laenge, j, T, generator, H, messdatei);//Schreibt Messwerte in Datei
	}
}

void messen(int laenge, double T, double j, int messungen, FILE *gitterdatei, FILE *messdatei, gsl_rng *generator){
	//Führt  messungen Messungen an Gitter in gitterdatei durch mit T, j, generator, speichert das Ergebnis in messdatei
	//generator für messen innerhalb der Funktion seeden?
	int gitter[laenge*laenge];
	einlesen(gitter, laenge, gitterdatei);
	double H=hamiltonian(gitter, laenge, j);
	for (int messung=0; messung<messungen; messung+=1){
		fprintf(messdatei,"%f\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		H=sweep(gitter, laenge, j, T, generator, H, messdatei);//Schreibt Messwerte in Datei
	}
	//printf("%d\n", gittersumme(gitter, laenge));
}

double mittelwertberechnung(FILE *messdatei, int messungen, const int spalte){
	//berechnet den Mittelwert aus spalte aus den messungen in messdatei
	double summe=0;//Speichert Summe über Messungen
	double einwert=0;//Speichert einen ausgelesenen Wert
	double ergebnisarray[3];//speichert alle messungen
	rewind(messdatei);//sichergehen, dass alle Messdaten verwendet werden
	for (int messung=0; messung<messungen; messung+=1){//Mittelwert über Messung bilden
		fscanf(messdatei, "%lf \t %lf \t %lf \n", &ergebnisarray[0], &ergebnisarray[1], &ergebnisarray[2]);
		einwert=ergebnisarray[spalte];//wählt korrekte messung aus
		summe+=einwert;
	}
	return summe/(double)messungen;
}

int main(int argc, char **argv){
	int laenge=300;
	double j=1.0;
	int seed=5;
	//double T=0.5;
	int messungen=1000;
	FILE *gitterthermdatei, *messdatei, *mittelwertdatei;
	const int temperaturzahl=15;
	double temperaturarray[15]={0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10};
	
	char dateinametherm[20], dateinamemessen[20];
	double mittelwert, varianz;
	
	gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	mittelwertdatei=fopen("messenmittel.txt", "w");//speichert mittlewerte ans Ende der Datei
	for (int n=0; n<temperaturzahl; n+=1){    //counting through given temperaturs
	 
		sprintf(dateinametherm,"thermalisierung-%.2d.txt",n);//.2, damit alle dateinamengleich lang sind
		sprintf(dateinamemessen,"messung-%.2d.txt",n);//.2, damit alle dateinamengleich lang sind
		gitterthermdatei = fopen(dateinametherm, "r");
		messdatei = fopen(dateinamemessen, "r");
		gsl_rng_set(generator, seed);//initialisieren, bei jedem Durchlauf mit gleichem seed
		//thermalisieren(laenge, temperaturarray[n], j, seed, gitterthermdatei, generator);
		//messen(laenge, temperaturarray[n], j, messungen, gitterthermdatei, messdatei, generator);
		mittelwert=mittelwertberechnung(messdatei, messungen);
		varianz=varianzberechnung(messdatei, messungen, mittelwert);
		fprintf(mittelwertdatei, "%f\t%f\t%f\n", temperaturarray[n], mittelwert, varianz);
		printf("%f\t%f\t%f\n", temperaturarray[n], mittelwert, varianz);
		
	}

	//~ gittertherm=fopen("gitterthermalisiert.txt", "w+");//speichert thermalisiertes Gitter, auch Auslesen möglich
	//~ thermalisieren(laenge, T, j, seed, gittertherm, generator);
	
	//~ //messen
	//~ messdatei=fopen("messdatei.txt", "w+");//speichert messungen, bei Bedarf auf w+ setzen!
	//~ messen(laenge, T, j, messungen, gittertherm, messdatei, generator);
	
	//~ mittelwertdatei=fopen("messenmittel.txt", "a");//speichert mittlewerte ans Ende der Datei
	//~ printf("%f\t%f\n", T, mittelwert(messdatei, messungen));	
	fclose(mittelwertdatei);
	fclose(messdatei);
	fclose(gitterthermdatei);
	
	gsl_rng_free(generator);//free, close: zum Verhindern von Speicherproblemen

	return 0;
}

//Metropolis:
//N_0 sweeps zum thermalisieren
//Wann thermalisiert? Schwelle für Änderung Hamiltonian?
//Plot: Änderung Hamiltonian wird nicht geringer für T=10
//Fuer T kleiner, kurz große Änderungen, aber danach imemr noch "grosse" Schwankungen
//Fuer groessere laengen Effekt sichtbarer, Abbruchkriterium: Diff=0
//N sweeps zum Messen: Was als Rückgabewert von sweep? sweep überhaupt als Funktion oder in main integrieren?
//Als erstes messen: Akzeptanzrate
//Wie Mittelwert über Messungen bilden?/Wie Messungen auswerten
//Ziel: Aktzeptanzrate(T) darstellen
//https://stackoverflow.com/questions/5922579/c-ignoring-a-comment-line-in-input-file
//https://cboard.cprogramming.com/c-programming/132421-need-help-opening-multiple-files-c-programming.html
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

	//double gamma_T=10.0; //benutzt für SA
	//int numberofsweeps, temperaturechanges;
//~ int gitter [laenge*laenge];//initialisiert gitter
	//~ initialisierung(gitter, laenge, seed);
	//~ //ausgabe(gitter, laenge);
	//~ double H=hamiltonian(gitter, laenge, j); //Anfangsenergie
	//~ printf("H=%f\n", H);
	
	//~ double Hneu=H;
	//~ double Halt=H+laenge*j+1;
	//~ //int sumdifferences;
	//~ //thermalisieren
	//~ FILE *therm=fopen("thermalisierung.txt", "w");
	//~ FILE *dummyfile=fopen("dummy.txt", "w");//speichert messergebnisse waehrend des thermalisierens->Nicht benötigt
	//~ FILE *gittertherm=fopen("gitterthermalisiert.txt", "w+");//speichert thermalisiertes Gitter, auch Auslesen möglich
	//~ int N0=0;//zählt, wie viele sweeps zum Thermalisieren benoetigt werden
	//~ while (Halt-Hneu>0){//Abbruchkriterium
		//~ Halt=Hneu;//Zustand der vorherigen Iteration speichern zum Vergleich
		//~ Hneu=sweep(gitter, laenge, j, T, generator, Halt, dummyfile);//neuen Zustand durch sweep vom alten Zustand
		//~ fprintf(therm,"%d %f\n", N0, Halt-Hneu);
		//~ N0+=1;
	//~ }
	//~ fclose(therm);
	//~ fclose(dummyfile);
	//~ ausgabe(gitter, laenge, gittertherm);//ermöglicht Ändern des verwendeten Gitter
//einlesen(gitter, laenge, gittertherm);
//fprintf(messen, "Messung\tVeraenderungen\tAkzeptanzrate\n");
	//~ for (int messung=0; messung<messungen; messung+=1){
		//~ fprintf(messdatei,"%d\t", messung);//Schreibt in Datei, um die wievielte Messung es sich handelt
		//~ H=sweep(gitter, laenge, j, T, generator, H, messdatei);//Schreibt Messwerte in Datei
	//~ }
	//~ double summe=0;//Speichert Summe über Messungen
	//~ double einwert=0;//Speichert einen ausgelesenen Wert
	//~ int dummywert;//zum Zuweisen nicht benoetigter Messdaten
	//~ rewind(messdatei);
	//~ for (int messung=0; messung<messungen; messung+=1){//Mittelwert über Messung bilden
		//~ fscanf(messdatei, "%d \t %d \t %lf \n", &dummywert, &dummywert, &einwert);
		//~ summe+=einwert;
	//~ }
