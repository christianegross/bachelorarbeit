//Christiane, start 15.04.20
//Erste Implementierung Ising-Modell für Bachelorarbeit
//Versuch Metropolis-Algorithmus


#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include <sys/time.h>//Zur Messung der Wallclocktime beim messen ->Vergleich der Sweep-Funktionen

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


double deltah(int *gitter, int d1, int d2, int laenge, double j){
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

void flipspin(int *gitter, int d1, int d2, int laenge){
	gitter[laenge*d1+d2]*=-1;
}
	
double sweepaltohnepar(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen){
	//geht das ganze Gitter durch und versucht, jeden Spin umzudrehen. Zählt die Veränderungen, misst Akzeptanzrate und Magneitiserung und gibt aktuellen Hamiltonian zurück
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
			delta=deltah(gitter, d1, d2, laenge, j);
			if (((d1+d2)%2==0)&&(tryflip(gitter, d1, d2, laenge, j, T, generator, delta)==1)){//Wenn schwarzer Punkt und Spin geflippt wurde
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
			delta=deltah(gitter, d1, d2, laenge, j);
			if (((d1+d2)%2==1)&&(tryflip(gitter, d1, d2, laenge, j, T, generator, delta)==1)){//Wenn weisser Punkt und Spin geflippt wurde
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

double sweep(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen){
	//geht erst alle schwarzen und dann alle weissen Punkte des Gitters durch, macht ein Metropolis-Update an jedem Punkt, schreibt Akzeptanzrate und MAgnetisierung in dateimessungen
	//arbeitet parallel in schleifen ueber die einzelnen Farben
	double H=hamiltonian;//misst Gesamtveraenderung
	double veraenderungH=0;//misst Veraenderung in einem parallen Thread
	double delta=0;
	int changes =0;//misst Gesamtzahl der spinflips
	int changesklein=0;//misst Spinflips in parallelen Thread
	//schwarz: d1+d2 gerade
	#pragma omp parallel firstprivate (delta, veraenderungH, changesklein) shared (H, changes)
	{
		#pragma omp for
		for (int d1=0; d1<laenge;d1+=1){
			for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				delta=deltah(gitter, d1, d2, laenge, j);
				if (((d1+d2)%2==0)&&(tryflip(gitter, d1, d2, laenge, j, T, generator, delta)==1)){//Wenn schwarzer Punkt und Spin geflippt wurde
					flipspin(gitter, d1, d2, laenge);//in Gitter speichern
					veraenderungH+=delta;//Zwischenvariable, damit es keine Konflikte beim updaten gibt
					changesklein+=1;
				}
			}
		}
		#pragma omp critical//damit das updaten keine konflikte verursacht
		{H+=veraenderungH;
			changes+=changesklein;}
	}
	veraenderungH=0;//Zuruecksetzen, damit in naechster paralleler Region nur deren Veraenderungen gezaehlt werden
	changesklein=0;
	//weiss: d1+d2 ungerade, sonst analog zu schwarz
	#pragma omp parallel firstprivate (delta, veraenderungH, changesklein) shared (H, changes)
	{
		#pragma omp for
		for (int d1=0; d1<laenge;d1+=1){
			for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
				delta=deltah(gitter, d1, d2, laenge, j);
				if (((d1+d2)%2==1)&&(tryflip(gitter, d1, d2, laenge, j, T, generator, delta)==1)){//Wenn weisser Punkt und Spin geflippt wurde
					flipspin(gitter, d1, d2, laenge);//in Gitter speichern
					veraenderungH+=delta;
					changesklein+=1;
				}
			}
		}
		#pragma omp critical
		{H+=veraenderungH;
		changes+=changesklein;}
	}
	double akzeptanzrate=(double)changes/(double)laenge/(double)laenge;
	double magnetisierung=(double)gittersumme(gitter, laenge)/(double)laenge/(double)laenge;
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



double mittelwertberechnungnaiv(FILE *messdatei, int messungen, const int spalte, const int spalten){
	//berechnet den Mittelwert aus spalte aus den messungen in messdatei, messdatei muss nur aus doubles in spalten bestehen
	double summe=0;//Speichert Summe über Messungen
	double einwert=0;//Speichert einen ausgelesenen Wert
	double ergebnisarray[spalten];//speichert alle messungen
	rewind(messdatei);//sichergehen, dass alle Messdaten verwendet werden
	for (int messung=0; messung<messungen; messung+=1){//Mittelwert über Messung bilden
		for (int i=0; i<spalten; i+=1){
			fscanf(messdatei, "%lf", &ergebnisarray[i]);//scannt einzelne doubles
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
			fscanf(messdatei, "%lf", &ergebnisarray[i]);
			if (i==spalten-1){fscanf(messdatei, "\n");}
		}
		einwert=ergebnisarray[spalte];//wählt korrekte messung aus
		summe+=(einwert-mittelwert)*(einwert-mittelwert);
	}
	return sqrt(summe/((double)messungen-1));
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
				fscanf(messdatei, "%lf", &ergebnisarray[i]);
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
	double *bootstraparray=(double*)malloc(sizeof(double)*r);//speichert ausgewählten Daten
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
	fprintf(ausgabedatei, "2\t%d\t%d\t%f\t%e\t%f\n", l,r, mittelwert, varianz, temperatur);//Ausgabe
	free(bootstraparray);
}

void bootstrap(int l, int r, int messungen, double temperatur, double *blockarray, gsl_rng *generator, FILE *ausgabedatei){
//berechnet Mittelwert und Varianz aus r gebootstrappten replikas, schreibt es in ausgabedatei
	double mittelwert=0;//speichern zwischenwerte
	double replica;//einzelen bootstrapreplica
	double varianz=0;
	double *bootstraparray=(double*)malloc(sizeof(double)*r);//speichert ausgewählten Daten
	#pragma omp parallel for private (replica)
	for (int durchgang=0; durchgang<r; durchgang+=1){//Zieht r replicas, speichert sie in bootstraparray
		replica=bootstrap_replication(l, messungen, blockarray, generator);
		bootstraparray[durchgang]=replica;//speichern fuer Varianzbildung
	}
	for (int durchgang=0; durchgang<r; durchgang+=1){//Mittelwert ueber gezogene Replikas
		mittelwert+=bootstraparray[durchgang];
	}
	mittelwert/=r;//Standardschaetzer
	for (int durchgang=0; durchgang<r; durchgang+=1){//Berechnet Varianz von gezogenen Replikas
		varianz+=(bootstraparray[durchgang]-mittelwert)*(bootstraparray[durchgang]-mittelwert);
	}
	varianz=sqrt(varianz/((double)r-1));//Standardschaetzer
	fprintf(ausgabedatei, "1.0\t%f\t%f\t%f\t%e\t%f\n", (double)l,(double)r, mittelwert, varianz, temperatur);//Ausgabe, 1, um zu zeigen, dass parallel gerechnet wurde
	free(bootstraparray);
}

void ableitung(int l, int temperaturen, const int spalten, const int spaltemessung, const int spaltefehler, const int spaltetemperatur, const int spaltel, FILE *messdatei, FILE *ausgabedatei){
	//Bildet die ableitung nach Zei-Punkt-Formel
	//Fehler nach Gaussscher Fehlerfortpflanzung
	double x1, x2, y1, y2, dy1, dy2;//Zur Berechnung benoetigte Groessen
	double mitteltemperatur, ableitung, fehlerableitung;//Groessen die berechnet werden sollen
	double ergebnisarray[spalten];//Speichert eingelesene Dateien
	rewind(messdatei);
	for (int i=0; i<spalten; i+=1){//erste Zeile scannen: Noch keine Ableitung moeglich
		fscanf(messdatei, "%lf", &ergebnisarray[i]);
		if (i==spalten-1){fscanf(messdatei, "\n");}
	}
	x1=ergebnisarray[spaltetemperatur];//Werte fuer erste Ableitung zuweisen
	y1=ergebnisarray[spaltemessung];	
	dy1=ergebnisarray[spaltefehler];
	for (int messung=1; messung<temperaturen; messung+=1){//Alle Zeilen durchgehen
		for (int i=0; i<spalten; i+=1){
			fscanf(messdatei, "%lf", &ergebnisarray[i]);
			if (i==spalten-1){fscanf(messdatei, "\n");}
		}
		if (ergebnisarray[spaltel]==l){//nur Ableitung berechnen, wenn l richtig ist
			x2=ergebnisarray[spaltetemperatur];//Wertezuweisen
			y2=ergebnisarray[spaltemessung];	
			dy2=ergebnisarray[spaltefehler];
			mitteltemperatur=(x1+x2)/2;//Berechnung der Temperatur, bei der die Ableitung berechnet wird
			ableitung=(y2-y1)/(x2-x1);//Zwei-Punkt-Formel mit variablem Abstand möglich
			fehlerableitung=(sqrt(dy1*dy1+dy2*dy2))/(x2-x1);//Gausssche Fehlerfortpflanzung
			fprintf(ausgabedatei, "%f\t%f\t%f\n", mitteltemperatur, ableitung, fehlerableitung);
			x1=x2;//Zuweisung fuer naechste Zeile
			y1=y2;
			dy1=dy2;
		}
	}	
}

int main(int argc, char **argv){
	//benoetigte Variablen initialisieren
	int anzahlcores=1;
	omp_set_num_threads(anzahlcores);//Setzt die nummer an Kernen, die in den parallelen Regionen verwendet werden.
	int laenge=50;//laenge der verwendeten Gitter
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int N01=100000;//sweeps beim ersten Thermalisieren
	int N0=10000;//benoetigte sweeps zum Thermalisieren
	int messungen=10000;//pro temperatur, zweierpotenz um blocken einfacher zu machen
	int r;//Anzahl an samples für den Bootstrap
	FILE *gitterthermdatei, *messdatei, *mittelwertdatei, *dummydatei, *bootstrapalledatei, *ableitungdatei, *zeitdatei;//benoetigte Dateien zur Ausgabe
	int temperaturzahl=300;//Temperaturen, beid enen gemessen wird
	char dateinametherm[150], dateinamemessen[150], dateinamemittel[150], dateinamebootstrapalle[150], dateinameableitung[150], dateinamezeit[150];//Um Dateien mit Variablen benennen zu koennen
	double mittelwertmag, varianzmag, mittelwertakz, varianzakz;//fuer naive Fehler
	double *temperaturarray=(double*)malloc(sizeof(double)*temperaturzahl);//speichert verwendete Temperaturen
	for (int i=0; i<temperaturzahl;i++){//Temperaturarray intalisieren
		temperaturarray[i]=0.015*i+0.015;
	}
	int l;//Laenge der Blocks
	double *blockarray;//Zum Speichern der geblockten Messwerte
	double blocklenarray[12]={32, 64,128, 256, 384, 512, 640, 758, 876, 1024, 1280, 1536};//Blocklaengen, bei denen gemessen wird
	gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	sprintf(dateinamemittel,"/hiskp4/gross/bachelorarbeit/Messungen/Mittelwerte/messenmittel-l%.4d-m-%.6d.txt",laenge, messungen);//speichert naive Mittelwerte
	sprintf(dateinamebootstrapalle,"/hiskp4/gross/bachelorarbeit/Messungen/Bootstrapges/bootstrapalle-l%.4d-m-%.6d.txt",laenge, messungen);//speichert Mitteelwerte aus Bootstrap
	sprintf(dateinameableitung,"/hiskp4/gross/bachelorarbeit/Messungen/ableitung-laenge-%.4d-m-%.6d.txt",laenge, messungen);//speichert Ableitung
	sprintf(dateinamezeit,"/hiskp4/gross/bachelorarbeit/Messungen/zeiten-laenge-%.4d-m-%.6d-cores-%.2d.txt",laenge, messungen, anzahlcores);//speichert Ableitung
	mittelwertdatei=fopen(dateinamemittel, "w+");
	bootstrapalledatei=fopen(dateinamebootstrapalle, "w+");
	ableitungdatei=fopen(dateinameableitung, "w");
	zeitdatei=fopen(dateinamezeit, "w");
	//Messen der zeit, die während des Programms vergeht, aus C-Kurs kopiert:
	struct timeval anfangmessen, endemessen, anfangbootstrap, endebootstrap, anfangprogramm, endeprogramm;
	double sec, usec, zeitgesmessen, summezeitgesmessen, zeitgesbootstrap, summezeitgesbootstrap, zeitgesprogramm;
	summezeitgesmessen=0;//Zeit fuer alle Temperaturen insgesamt
	summezeitgesbootstrap=0;//Zeit fuer alle Temperaturen insgesamt
	gettimeofday(&anfangprogramm, NULL);
	//Thermalisierung des ersten Gitters, nicht ueber letztes verwendetes Gitter moeglich
	int gitter[laenge*laenge];
	initialisierung(gitter, laenge, seed);
	dummydatei=fopen("dummytherm.txt", "w");//speichert Gitter nach dem ersten Thermalisieren, das nicht benutzt wird
	gsl_rng_set(generator, seed);
	thermalisieren(laenge, temperaturarray[0], j, seed, N01, gitter, dummydatei, generator);//Erstes Thermalisierens, Anzahl je nach Länge groesser machen
	fclose(dummydatei);
	for (int n=0; n<temperaturzahl; n+=10){    //ueber alle gegebenen Temperaturen messen
		//printf("%d\n", n);
		sprintf(dateinametherm,"/hiskp4/gross/bachelorarbeit/Messungen/ThermalisierteGitter/thermalisierung-laenge%.4d-m%.6d-t%.3d.txt",laenge,messungen,n);//.2, damit alle dateinamengleich lang sind
		sprintf(dateinamemessen,"/hiskp4/gross/bachelorarbeit/Messungen/Messwerte/messung-laenge%.4d-m%.6d-t%.3d.txt",laenge,messungen,n);//.2, damit alle dateinamengleich lang sind
		gitterthermdatei = fopen(dateinametherm, "w+");//Zum speichern der thermalisierten Gitter
		messdatei = fopen(dateinamemessen, "w+");//Zum Speichern der Messdaten
		gsl_rng_set(generator, seed);//initialisieren, bei jedem Durchlauf mit gleichem seed
		thermalisieren(laenge, temperaturarray[n], j, seed, N0, gitter, gitterthermdatei, generator);
		gettimeofday(&anfangmessen, NULL);
		messen(laenge, temperaturarray[n], j, messungen, gitterthermdatei, messdatei, generator);
		gettimeofday(&endemessen, NULL);
		sec= (double)(endemessen.tv_sec-anfangmessen.tv_sec);
		usec= (double)(endemessen.tv_usec-anfangmessen.tv_usec);
		zeitgesmessen=sec+1e-06*usec;
		summezeitgesmessen+=zeitgesmessen;
		//printf("bei T=%f haben %d Messungen %f Sekunden gebraucht\n", temperaturarray[n], messungen, zeitgesmessen);
		fprintf(zeitdatei, "0.0\t%f\t%f\t%f\n", temperaturarray[n], (double)messungen, zeitgesmessen);
		mittelwertakz=mittelwertberechnungnaiv(messdatei, messungen, 1, 3);
		varianzakz=varianzberechnungnaiv(messdatei, messungen, mittelwertakz, 1, 3);
		mittelwertmag=mittelwertberechnungnaiv(messdatei, messungen, 2, 3);
		varianzmag=varianzberechnungnaiv(messdatei, messungen, mittelwertmag, 2, 3);
		fprintf(mittelwertdatei, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", (double)laenge, temperaturarray[n],j/temperaturarray[n], mittelwertakz, varianzakz, mittelwertmag, varianzmag);
		gettimeofday(&anfangbootstrap, NULL);
		for(int len=0;len<12;len+=1){//Fuer verschiedene l blocking und bootstrapping durchfuehren
			l=blocklenarray[len];
			//printf("%d\t%d\n", n, l);
			blockarray=(double*)malloc(sizeof(double)*messungen/l);//zum Speichern der Blocks
			r=4*messungen;//Anzahl an Replikas, die beim Bootstrappen erzeugt werden
			blocks_generieren(l, messungen, 2, 3, blockarray, messdatei);//blocking
			//Vergleich bootstrapping mit und ohne parallelisierung
			bootstrap(l, r, messungen, temperaturarray[n], blockarray, generator,bootstrapalledatei);//bootstrapping
			//bootstrapohnepar(l, r, messungen, temperaturarray[n], blockarray, generator,bootstrapalledatei);//bootstrapping
			free(blockarray);
		}
		gettimeofday(&endebootstrap, NULL);
		sec= (double)(endebootstrap.tv_sec-anfangbootstrap.tv_sec);
		usec= (double)(endebootstrap.tv_usec-anfangbootstrap.tv_usec);
		zeitgesbootstrap=sec+1e-06*usec;
		summezeitgesbootstrap+=zeitgesbootstrap;
		fprintf(zeitdatei, "1.0\t%f\t%f\t%f\n", temperaturarray[n], (double)messungen, zeitgesbootstrap);
		//printf("bei T=%f hat das Bootstrapping %f Sekunden gebraucht\n", temperaturarray[n], zeitgesbootstrap); 
		fclose(messdatei);
		fclose(gitterthermdatei);
	}
	ableitung(128, 30*12, 6, 3,4,5,1, bootstrapalledatei, ableitungdatei);
	gettimeofday(&endeprogramm, NULL);
	sec= (double)(endeprogramm.tv_sec-anfangprogramm.tv_sec);
	usec= (double)(endeprogramm.tv_usec-anfangprogramm.tv_usec);
	zeitgesprogramm=sec+1e-06*usec;
	printf("Insgesamt hat das Messen %f Sekunden gebraucht\n", summezeitgesmessen);
	printf("Insgesamt hat das Botstrapping %f Sekunden gebraucht\n", summezeitgesbootstrap);
	printf("Insgesamt hat das Programm %f Sekunden gebraucht\n", zeitgesprogramm);
	fclose(mittelwertdatei);
	fclose(bootstrapalledatei);
	fclose(ableitungdatei);
	fclose(zeitdatei);
	free(temperaturarray);
	gsl_rng_free(generator);//free, close: zum Verhindern von Speicherproblemen

	return 0;
}
