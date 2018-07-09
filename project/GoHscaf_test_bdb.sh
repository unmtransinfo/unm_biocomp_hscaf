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
#
java -classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds \
	-i $HOME/data/smi/hscaf_testset.smi \
	-o $HOME/Download/z_hscaf_out_nodb.smi \
	-out_scaf $HOME/Download/z_hscaf_scaf_nodb.smi \
	-inc_mol \
	-inc_scaf \
	-inc_link \
	-inc_chain \
	-scaflist_append2title \
	-nodb \
#
#
java -classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds \
	-i $HOME/data/smi/hscaf_testset.smi \
	-o $HOME/Download/z_hscaf_out_db.smi \
	-out_scaf $HOME/Download/z_hscaf_scaf_db.smi \
	-inc_mol \
	-inc_scaf \
	-inc_link \
	-inc_chain \
	-scaflist_append2title \
	-destroyexistingdb \
	-dbdir /tmp/hscaf \
	-vv
#
#
