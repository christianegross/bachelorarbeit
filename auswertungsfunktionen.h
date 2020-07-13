//Christiane, 25.05.20
//Header-Datei mit den Funktionen, die zum Auswerten in der Bachelorarbeit benoetigt werden
#pragma once


double mittelwertberechnungnaiv(FILE *messdatei, int messungen, const int spalte, const int spalten);
double varianzberechnungnaiv(FILE *messdatei, int messungen, double mittelwert, const int spalte, const int spalten);
void blocks_generieren(int l, int messungen, const int spalte, const int spalten, double *blockarray,  FILE *messdatei);
double bootstrap_replication(int l, int messungen, double *blockarray, gsl_rng *generator);
void bootstrapohnepar(int l, int r, int messungen, double temperatur, double *blockarray, gsl_rng *generator, FILE *ausgabedatei);
void bootstrap(int l, int r, int messungen, double temperatur, double *blockarray, gsl_rng **generatoren, FILE *ausgabedatei);
void ableitung(int l, int temperaturen, const int spalten, const int spaltemessung, const int spaltefehler, const int spaltetemperatur, const int spaltel, FILE *messdatei, FILE *ausgabedatei);
void ableitungdreipunkt(int l, int temperaturen, const int spalten, const int spaltemessung, const int spaltefehler, const int spaltetemperatur, const int spaltel, FILE *messdatei, FILE *ausgabedatei);
double mittelwertarray(double *array, int messungen);
double varianzarray(double *array, int messungen, double mittelwert);
double minarray(double *array, int messungen);void minline(FILE *datei, const int spalten, const int minwertspalte, const int zeilen);
void sqrtspalte(FILE *einlesedatei, FILE *ausgabedatei, const int spalten, const int spalte, const int zeilen);
void maxline(FILE *datei, const int spalten, const int maxwertspalte, const int zeilen);
