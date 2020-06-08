//Christiane, 25.05.20
//Header-Datei mit den Funktionen, die zum Messen in der Bachelorarbeit benoetigt werden
#pragma once
void initialisierung(char *gitter, int laenge, int seed);
void ausgabemitplot(char *gitter, int laenge, FILE *datei, FILE *plotdatei);
void ausgabe(char *gitter, int laenge, FILE *datei);
void einlesen(char *gitter, int laenge, FILE *datei);
int gittersummeohnepar (char *gitter, int laenge);
int gittersumme (char *gitter, int laenge);
double hamiltonian(char *gitter, int laenge, double j);
double deltahalt(char *gitter, int d1, int d2, int laenge, double j);
int tryflipalt(char *gitter,  int d1, int d2, int laenge, double j, double T, gsl_rng *generator, double delta);
int deltah(char *gitter, int d1, int d2, int laenge);
int tryflip(char *gitter,  int d1, int d2, int laenge, double j, double T, gsl_rng *generator, double wahrscheinlichkeit);
void flipspin(char *gitter, int d1, int d2, int laenge);
double wahrscheinlichkeit(int delta, double *wahrscheinlichkeiten);
double sweepaltohnepar(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweepalt(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweep(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
void thermalisierenmitplot(int laenge, double T, double j, int seed,int N0, char *gitter, FILE *ausgabedatei, FILE *plotdatei, gsl_rng *generator);
void thermalisieren(int laenge, double T, double j, int seed,int N0, char *gitter, FILE *ausgabedatei, gsl_rng *generator);
void messen(int laenge, double T, double j, int messungen, FILE *gitterdatei, FILE *messdatei, gsl_rng *generator);
void messenvergleichen(int laenge, double T, double j, int messungen, FILE *gitterdatei, FILE *messdatei, FILE *vergleichsdatei, gsl_rng *generator);

