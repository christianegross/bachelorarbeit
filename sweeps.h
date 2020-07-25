//Christiane, 16.06.20
//header-Datei fuer die verschiedenen Sweep-Funktionen
//Erklaerungen in sweeps.c, nur sweepaltohnepar, sweepalt und sweepmehreregneratoren benutzt, der Rest sind Zwischenversionen
#pragma once


double sweepaltohnepar(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweepalt(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweepmehreregeneratoren(char *gitter, int laenge, double j, double T, gsl_rng **generatoren, double hamiltonian, double *wahrscheinlichkeiten, FILE *dateimessungen);
double sweepzweipar(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweep(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweeplookup(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen,int *lookupplus, int *lookupminus);
double sweepmehreregeneratorenv0(char *gitter, int laenge, double j, double T, gsl_rng **generatoren, double hamiltonian, FILE *dateimessungen);
double sweepmehreregeneratorenv1(char *gitter, int laenge, double j, double T, gsl_rng **generatoren, double hamiltonian, FILE *dateimessungen);
double sweepmehreregeneratorenv2(char *gitter, int laenge, double j, double T, gsl_rng **generatoren, double hamiltonian, FILE *dateimessungen);
double sweepmehreregeneratorenv3(char *gitter, int laenge, double j, double T, gsl_rng **generatoren, double hamiltonian, FILE *dateimessungen);
double sweepmehreregeneratorenv4(char *gitter, int laenge, double j, double T, gsl_rng **generatoren, double hamiltonian, double *wahrscheinlichkeiten, FILE *dateimessungen);
double sweepeineschleife(char *gitter, int laenge, double j, double T, gsl_rng **generatoren, double hamiltonian, FILE *dateimessungen);

