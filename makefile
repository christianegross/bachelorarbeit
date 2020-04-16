ising: ising.o
	gcc -std=c99 -Wall -pedantic -o ising ising.o -lgsl -lgslcblas -lm
	
ising.o: ising.c
	gcc -std=c99 -Wall -pedantic -I /usr/include/ ising.c -c
