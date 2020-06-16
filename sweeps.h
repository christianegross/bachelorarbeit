//Christiane, 16.06.20
//header-Datei fuer die verschiedenen Sweep-Funktionen
#pragma once


double sweepaltohnepar(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweepalt(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweepzweipar(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweep(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen);
double sweepmehreregeneratoren(char *gitter, int laenge, double j, double T, gsl_rng **generatoren, double hamiltonian, FILE *dateimessungen);
double sweeplookup(char *gitter, int laenge, double j, double T, gsl_rng *generator, double hamiltonian, FILE *dateimessungen,int *lookupplus, int *lookupminus);

