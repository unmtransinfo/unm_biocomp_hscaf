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
time java -classpath $CLASSPATH:$HOME/src/java/lib/unm_biocomp_hscaf.jar \
	edu.unm.health.biocomp.hscaf.hier_scaffolds \
	-i $HOME/chemaxon/hscaf/data/huge_hiers_trees.smi \
	-out_scaf $HOME/chemaxon/hscaf/data/huge_hiers_trees_scaf.smi \
	-vv
#
#
time java -classpath $CLASSPATH:/home/app/lib/unm_biocomp_hscaf.jar \
	edu.unm.health.biocomp.hscaf.hier_scaffolds \
	-i $HOME/chemaxon/hscaf/data/huge_hiers_trees.smi \
	-out_scaf $HOME/chemaxon/hscaf/data/huge_hiers_trees_scaf.smi \
	-vv
#
#
