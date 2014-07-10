#!/bin/sh

BASEDIR=CHANGE.THIS

mkdir $BASEDIR/infiles 
mkdir $BASEDIR/outfiles
mkdir $BASEDIR/PoDA-out
rm $BASEDIR/nashPoDA.swarm

for i in `seq 0 100`
	do for p in `seq 1 36`
		do echo "pathBlock=${p}; i=${i}; setwd(\"$BASEDIR\"); source(\"$BASEDIR/nashPoDA.skel.R\")" > $BASEDIR/infiles/nashPoDA.pathBlock_${p}.run_${i}.R
		echo "/usr/local/bin/R --vanilla < $BASEDIR/infiles/nashPoDA.pathBlock_${p}.run_${i}.R > $BASEDIR/outfiles/nashPoDA.pathBlock_${p}.run_${i}.Rout" >> $BASEDIR/nashPoDA.swarm
	done
done

echo "submit with: swarm -f nashPoDA.swarm"

