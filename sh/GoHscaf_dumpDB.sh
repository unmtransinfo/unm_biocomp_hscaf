#!/bin/sh
#
set -x
#
if [ `uname -s` = "Darwin" ]; then
	appdir=/Users/app
else
	appdir=/home/app
fi
#
CLASSPATH=$HOME/src/java/lib/unm_biocomp_hscaf.jar
#CLASSPATH=$appdir/lib/unm_biocomp_hscaf.jar
CLASSPATH=$CLASSPATH:$appdir/lib/jchem.jar
CLASSPATH=$CLASSPATH:$appdir/lib/berkeleydb.jar
#
HOST=`hostname -s`
#
scratchdir=$HOME/projects/badapple/scratch.$HOST
#scratchdir=scratch
#
java -classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds \
	-dumpdb data/dumpdb_scaf.smi \
	-dbdir $scratchdir \
	-vv
#
wc -l data/dumpdb_scaf.smi
#
