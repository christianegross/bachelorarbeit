//Christiane, 10.07.20
//Header-Datei mit Funktionen, die fuer Programm mit MPI benoetigt werden
//Viele Funktionen aus entsprechenden OpenMP-Codes kopiert
//Erklaerungen zu Funktionen in mpifunktionen.c
#pragma once
int deltahmpi(int d1, int d2, int laenge, int teillaenge, char *untergitter, char *nachbarunten, char *nachbaroben);
void initialisierenhomogen(char *gitter, int laenge);
void einlesen(char *gitter, int laenge, FILE *datei);
void ausgabe(char *gitter, int laenge, FILE *datei);
double hamiltonian(char *gitter, int laenge, double j);
int tryflip(gsl_rng *generator, double wahrscheinlichkeit);
double mittelwertberechnungnaiv(FILE *messdatei, int messungen, const int spalte, const int spalten);
double varianzberechnungnaiv(FILE *messdatei, int messungen, double mittelwert, const int spalte, const int spalten);
double sweepmpi(int laenge, FILE *ausgabedatei, gsl_rng **generatoren, int anzproz, int myrank, int teillaenge, char *untergitter, char *nachbarunten, char *nachbaroben, double* wahrscheinlichkeiten, double j, double H);
void messenmpi(int messungen, int laenge, double T, double j, char *gitter, FILE *ausgabedatei, gsl_rng **generatoren);
void thermalisierenmpi(int messungen, int laenge, double T, double j, char *gitter, FILE *dummymessung, FILE *ausgabedatei, gsl_rng **generatoren);
void blocks_generieren(int l, int messungen, const int spalte, const int spalten, double *blockarray,  FILE *messdatei);
double bootstrap_replication(int l, int messungen, double *blockarray, gsl_rng *generator);
void bootstrapohnepar(int l, int r, int messungen, double temperatur, double *blockarray, gsl_rng *generator, FILE *ausgabedatei);
double mittelwertarray(double *array, int messungen);
double varianzarray(double *array, int messungen, double mittelwert);
double minarray(double *array, int messungen);void minline(FILE *datei, const int spalten, const int minwertspalte, const int zeilen);
double maxarray(double *array, int messungen);
