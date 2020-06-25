//Christiane, ab 25.05.20
//Funktionen fuer die Bachelorarbeit, die zum Auswerten der Messdaten benoetigt werden

#include <stdio.h>
#include <gsl/gsl_rng.h>//Zufallszahlen
#include "math.h"//exp-Funktion
#include <omp.h>//Parallelisierung
#include "auswertungsfunktionen.h"

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
	fprintf(ausgabedatei, "2\t%4d\t%d\t%f\t%e\t%f\n", l,r, mittelwert, varianz, temperatur);//Ausgabe
	free(bootstraparray);
}

void bootstrap(int l, int r, int messungen, double temperatur, double *blockarray, gsl_rng **generatoren, FILE *ausgabedatei){
//berechnet Mittelwert und Varianz aus r gebootstrappten replikas, schreibt es in ausgabedatei
	double mittelwert=0;//speichern zwischenwerte
	double replica;//einzelen bootstrapreplica
	double varianz=0;
	double *bootstraparray;
	if((bootstraparray=(double*)malloc(sizeof(double)*r))==NULL){//speichert ausgewählten Daten, prüft, ob Speicherplatz richitg bereitgestellt wurde
		printf("Fehler beim Allokieren des arrays fuer die Replika!\n");
		//return (-1);
	}
	#pragma omp parallel private (replica)
	{
	int threadnummer=omp_get_thread_num();
	#pragma omp for
	for (int durchgang=0; durchgang<r; durchgang+=1){//Zieht r replicas, speichert sie in bootstraparray
		replica=bootstrap_replication(l, messungen, blockarray, generatoren[threadnummer]);
		bootstraparray[durchgang]=replica;//speichern fuer Varianzbildung
	}
	}
	//~ mittelwert=mittelwertarray(bootstraparray, r);
	//~ varianz=varianzarray(bootstraparray, r, mittelwert);
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
