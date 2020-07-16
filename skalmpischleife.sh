#!/bin/bash

#for laenge in 12 24 36 48 72 96 120 144 192 384 768
#do
#	sbatch batchskriptmpiskalierungn6 ${laenge} 0.5 mpiskalt05n6
#	echo ${laenge}
#done

for temperatur in 1 2 3 4 42
do
	for anzahl in 1 2 3 4 6
	do
		sbatch batchskriptmpiskalierungn${anzahl} 768 temperatur mpiskalt${temperatur}
		echo ${temperatur}
		echo ${anzahl}
	done
done
