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
java -classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds_common \
	-q $HOME/chemaxon/hscaf/data/quinoline.sdf \
	-db $HOME/chemaxon/hscaf/data/nci_quinolines2d_input.sdf \
	-o $HOME/chemaxon/hscaf/data/common_scaf.smi \
	-vv
#
#
