#pragma once
void initialisierung(int *gitter, int laenge, int seed);
void ausgabe(int *gitter, int laenge, FILE *datei);
void einlesen(int *gitter, int laenge, FILE *datei);
int gittersummeohnepar (int *gitter, int laenge);
int gittersumme (int *gitter, int laenge);
double hamiltonian(int *gitter, int laenge, double j);
double deltahalt(int *gitter, int d1, int d2, int laenge, double j);
int tryflipalt(int *gitter,  int d1, int d2, int laenge, double j, double T, gsl_rng *generator, double delta);
int deltah(int *gitter, int d1, int d2, int laenge);
int tryflip(int *gitter,  int d1, int d2, int laenge, double j, double T, gsl_rng *generator, double wahrscheinlichkeit);
void flipspin(int *gitter, int d1, int d2, int laenge);
double wahrscheinlichkeit(int delta, double *wahrscheinlichkeiten);
double sweepaltohnepar(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweepalt(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweep(int *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
void thermalisieren(int laenge, double T, double j, int seed,int N0, int *gitter, FILE *ausgabedatei, gsl_rng *generator);
void messen(int laenge, double T, double j, int messungen, FILE *gitterdatei, FILE *messdatei, gsl_rng *generator);

