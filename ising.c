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

int gittersumme (int *gitter, int laenge){
	//berechnet Summe aller Elemente eines Gitter mit laenge*laenge
	int summe=0;
	for (int d1=0; d1<laenge; d1+=1){//geht in erster dimension durch (Zeile
		for (int d2=0; d2<laenge; d2+=1){//geht in zweiter dimension durch (alle Spalten einer Zeile)
			summe+=gitter[laenge*d1+d2];
		}
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
	fprintf(dateimessungen, "%f\t%f\n",akzeptanzrate, magnetisierung );//benoetigte messungen: Anzahl Veränderungen+Akzeptanzrate=Veränderungen/Möglichkeiten+Magnetisierung
	return H;
}

void thermalisieren(int laenge, double T, double j, int seed,int N0, int *gitter, FILE *ausgabedatei, gsl_rng *generator){
	//erzeugt ein thermalisiertes Gitter mit laenge*laenge, T, j, seed in ausgabedatei
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
	//generator für messen innerhalb der Funktion seeden?
	int gitter[laenge*laenge];
	einlesen(gitter, laenge, gitterdatei);
	double H=hamiltonian(gitter, laenge, j);
	for (int messung=0; messung<messungen; messung+=1){
		fprintf(messdatei,"%f\t", (double)messung);//Schreibt in Datei, um die wievielte Messung es sich handelt, double, damit Mittelwertbestimmung einfacher wird
		H=sweep(gitter, laenge, j, T, generator, H, messdatei);//Geht Gitter durch und schreibt Messwerte in Datei
	}
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


double varianzberechnung(FILE *messdatei, int messungen, double mittelwert, const int spalte){
	//berechnet die varianz über die gegebenen Messungen in messdatei mit mittelwert mittels Standardschaetzer
	double summe=0;//Speichert Summe über Messungen
	double einwert=0;//Speichert einen ausgelesenen Wert
	double ergebnisarray[3];//speichert alle messungen
	rewind(messdatei);//sichergehen, dass alle Messdaten verwendet werden
	for (int messung=0; messung<messungen; messung+=1){//Mittelwert über Messung bilden
		fscanf(messdatei, "%lf \t %lf \t %lf \n", &ergebnisarray[0], &ergebnisarray[1], &ergebnisarray[2]);
		einwert=ergebnisarray[spalte];//wählt korrekte messung aus
		summe+=(einwert-mittelwert)*(einwert-mittelwert);
	}
	return sqrt(summe/((double)messungen-1));
}

void blocks_generieren(int l, int messungen, const int spalte, double *blockarray,  FILE *messdatei, FILE *ausgabedatei){
	//Blockt die Daten aus spalte in Messdatei in Blöcke der Laenge l, Ergebnis in blockarray und ausgabedatei
	double zwischensumme, einwert;//speichern der Summe über die Messwerte und der einzelnen Werte
	double ergebnisarray[3];//speichern der ganzen Zeile aus der Messdatei
	rewind(messdatei);
	for (int block=0; block<messungen/l; block+=1){//jedes einzelne Element in Blockarray durchgehen
		zwischensumme=0;
		for (int wert=0; wert<l; wert+=1){//generiert einzelnes Element des blocks
			fscanf(messdatei, "%lf \t %lf \t %lf \n", &ergebnisarray[0], &ergebnisarray[1], &ergebnisarray[2]);
			einwert=ergebnisarray[spalte];//wählt korrekte messung aus
			zwischensumme+=einwert;
		}
		blockarray[block]=zwischensumme/(double)l;
		fprintf(ausgabedatei, "%f\n", blockarray[block]);
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

void bootstrap(int l, int r, int messungen, double temperatur, double *blockarray, gsl_rng *generator, FILE *ausgabedatei){
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
	fprintf(ausgabedatei, "%d\t%d\t%f\t%e\t%f\n", l,r, mittelwert, varianz, temperatur);//Ausgabe
	free(bootstraparray);
}

int main(int argc, char **argv){
	//benoetigte Variablen initialisieren
	int laenge=50;//laenge der verwendeten Gitter
	double j=1.0;
	int seed=5;//fuer den zufallsgenerator
	int messungen=4096*4;//pro temperatur, zweierpotenz um blocken einfacher zu machen
	int r;//Anzahl an samples für den Bootstrap
	FILE *gitterthermdatei, *messdatei, *mittelwertdatei, *dummydatei, *bootstrapdatei, *blockdatei, *bootstrapalledatei;//benoetigte Dateien zur Ausgabe
	int temperaturzahl=300;//Temperaturen, beid enen gemessen wird
	int N0=5000;//benoetigte sweeps zum Thermalisieren
	char dateinametherm[100], dateinamemessen[100], dateinamemittel[100], dateinamebootstrap[100], blockdateiname[100], dateinamebootstrapalle[100];//Um Dateien mit Variablen benennen zu koennen
	double mittelwertmag, varianzmag, mittelwertakz, varianzakz;
	double *temperaturarray=(double*)malloc(sizeof(double)*temperaturzahl);
	for (int i=0; i<temperaturzahl;i++){//Temperaturarray intalisieren
		temperaturarray[i]=0.015*i+0.015;
	}
	int l;//Laenge der Blocks
	double *blockarray;//Zum Speichern der geblockten Messwerte
	double blocklenarray[12]={32, 64,128, 256, 384, 512, 640, 758, 876, 1024, 1280, 1536};//Blocklaengen, bei denen gemessen wird
	gsl_rng *generator=gsl_rng_alloc(gsl_rng_mt19937);//Mersenne-Twister
	sprintf(dateinamemittel,"Messungen/Mittelwerte/messenmittel-l%.4d-m-%.6d.txt",laenge, messungen);//.2, damit alle dateinamengleich lang sind
	sprintf(dateinamebootstrapalle,"Messungen/bootstrapalle-l%.4d-m-%.6d.txt",laenge, messungen);//.2, damit alle dateinamengleich lang sind
	mittelwertdatei=fopen(dateinamemittel, "w");
	bootstrapalledatei=fopen(dateinamebootstrapalle, "w");
	//Thermalisierung des ersten Gitters, nicht ueber letztes verwendetes Gitter moeglich
	int gitter[laenge*laenge];
	initialisierung(gitter, laenge, seed);
	dummydatei=fopen("dummytherm.txt", "w");
	thermalisieren(laenge, temperaturarray[0], j, seed, 10000, gitter, dummydatei, generator);
	fclose(dummydatei);
	for (int n=0; n<temperaturzahl; n+=10){    //ueber alle gegebenen Temperaturen messen
		printf("%d\n", n);
		sprintf(dateinametherm,"Messungen/ThermalisierteGitter/thermalisierung-laenge%.4d-t%.3d.txt",laenge,n);//.2, damit alle dateinamengleich lang sind
		sprintf(dateinamemessen,"Messungen/Messwerte/messung-laenge%.4d-t%.3d.txt",laenge,n);//.2, damit alle dateinamengleich lang sind
		sprintf(dateinamebootstrap,"Messungen/Bootstrapwerte/bootstrap-laenge%.4d-t%.3d.txt",laenge,n);//.2, damit alle dateinamengleich lang sind
		
		gitterthermdatei = fopen(dateinametherm, "r+");//Zum speichern der thermalisierten Gitter
		messdatei = fopen(dateinamemessen, "r+");//Zum Speichern der Messdaten
		bootstrapdatei=fopen(dateinamebootstrap, "r+");//Zum Speichern der Werte, die beim Bootstrapping berechnet werden
		gsl_rng_set(generator, seed);//initialisieren, bei jedem Durchlauf mit gleichem seed
		//thermalisieren(laenge, temperaturarray[n], j, seed, N0, gitter, gitterthermdatei, generator);
		//messen(laenge, temperaturarray[n], j, messungen, gitterthermdatei, messdatei, generator);
		mittelwertakz=mittelwertberechnung(messdatei, messungen, 1);
		varianzakz=varianzberechnung(messdatei, messungen, mittelwertakz, 1);
		mittelwertmag=mittelwertberechnung(messdatei, messungen, 2);
		varianzmag=varianzberechnung(messdatei, messungen, mittelwertmag, 2);
		fprintf(mittelwertdatei, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", laenge, temperaturarray[n],j/temperaturarray[n], mittelwertakz, varianzakz, mittelwertmag, varianzmag);
		//fprintf(mittelwertdatei, "%d\t%f\t%f\n",n, mittelwertmag, varianzmag);
		//fprintf(bootstrapdatei, "l\tr\tmittelwert\tvarianz\n");
		for(int len=0;len<12;len+=1){//Fuer verschiedene l blocking und bootstrapping durchfuehren
			l=blocklenarray[len];
			//printf("%d\t%d\n", n, l);
			sprintf(blockdateiname,"Bootstraptest/Blocks/blocks-laenge%.4d-t%.3d-l%.6d.txt",laenge,n, l);//.2, damit alle dateinamengleich lang sind
			blockdatei=fopen(blockdateiname, "w+");
			blockarray=(double*)malloc(sizeof(double)*messungen/l);//zum Speichern der Blocks
			r=4*messungen;//Anzahl an Replikas, die beim Bootstrappen erzeugt werden
			blocks_generieren(l, messungen, 2, blockarray, messdatei, blockdatei);//blocking
			bootstrap(l, r, messungen, temperaturarray[n], blockarray, generator,bootstrapalledatei);//bootstrapping
			free(blockarray);
			fclose(blockdatei);
		}
		//printf("set title \"T=%f\"\n\nset ylabel \"Mittelwert\"\nf(x)=%f\n", temperaturarray[n], mittelwertmag);
		//printf("plot \"Messungen/Bootstrapwerte/bootstrap-laenge%.4d-t%.3d.txt\" u 1:3:4 w yerrorbars lt 7 ps 0.3 title \"Bootstrap-Mittelwerte\", f(x) lt 6 title \"naiver Mittelwert\"\n\n", laenge, n);
		//printf("set ylabel \"Varianz\"\nplot \"Messungen/Bootstrapwerte/bootstrap-laenge%.4d-t%.3d.txt\" u 1:4 lt 7 ps 0.4 title \"Bootstrap-Varianz\", \"Messungen/Bootstrapwerte/bootstrap-laenge%.4d-t%.3d.txt\" u 1:4 w lines lt 7 title \"Bootstrap-Varianz\"\n\n", laenge, n, laenge, n);
		fclose(messdatei);
		fclose(gitterthermdatei);
		fclose(bootstrapdatei);
	}
		
	fclose(mittelwertdatei);
	fclose(bootstrapalledatei);
	free(temperaturarray);
	gsl_rng_free(generator);//free, close: zum Verhindern von Speicherproblemen

	return 0;
}
//Zu klaeren: Wie oft generator neu initialisieren? Bei jeder Verwendung, damit die Ergebnisse reproduzierbar sind?

//Ergebnisse zum Mittelwert bilden in array funktioniert nicht, da unterschiedliche types in datei stehen ->einzelne Funktionen oder Änderung von sweep je nach Messung?
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
//funktion "thermalisieren" benötigt?
//Für große Gitter: thermalisieren dauert lange, Ergebnis in Datei speichern, nachher einlesen
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


//~ , *thermenergiedatei
//~ , dateinameenergie[20]
//~ sprintf(dateinameenergie,"thermenergie-%.2d.txt",n);//.2, damit alle dateinamengleich lang sind
//~ thermenergiedatei = fopen(dateinameenergie, "w+");
//~ thermalisierenenergie(laenge, temperaturarray[n], j, seed, thermenergiedatei, generator);//Zeigt Energieänderungen während des Thermalisierens
//~ printf("set title \"T=%.2f\" font \"arial,40\"\nplot \"thermenergie-%.2d.txt\" using 1:2 lt 7 ps 0.5 w lines\n\n", temperaturarray[n], n);

