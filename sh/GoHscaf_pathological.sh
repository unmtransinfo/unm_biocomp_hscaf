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
IFILE=$HOME/chemaxon/hscaf/data/hscaf_pathological.smi
PREFIX=$HOME/chemaxon/hscaf/data/z_pathological
#
java -classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds \
	-i $IFILE \
	-o ${PREFIX}_hscaf_out.smi \
	-out_scaf ${PREFIX}_hscaf_scaf.smi \
	-inc_mol \
	-inc_scaf \
	-inc_link \
	-inc_chain \
	-scaflist_append2title \
	-nodb \
	-maxrings 5 \
	-vv
#
#	-dbdir scratch -destroyexistingdb \
###
