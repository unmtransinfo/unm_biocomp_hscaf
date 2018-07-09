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
CLASSPATH=$CLASSPATH:/home/app/lib/unm_biocomp_util.jar
CLASSPATH=$CLASSPATH:$appdir/lib/jchem.jar
CLASSPATH=$CLASSPATH:$appdir/lib/berkeleydb.jar
#
#INFILE=$HOME/Download/mlsmr_100.smi
#INFILE=$HOME/Download/mlsmr_5000.smi
INFILE=$HOME/Download/mlsmr_1000.smi
#INFILE=data/GL0021.smi
#INFILE=data/testset.smi
#
#time java -classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds \
#	-i $INFILE \
#	-o data/mlsmr_hscaf_out.smi \
#	-out_scaf data/mlsmr_hscaf_scaf.smi \
#	-inc_mol -inc_scaf
#
CLASSPATH=$CLASSPATH_ORIG:$HOME/src/java/lib/unm_biocomp_hscaf.jar:$HOME/src/java/lib/unm_biocomp_util.jar
#
time java -classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds \
	-i $INFILE \
	-o data/mlsmr_hscaf_out_dev.smi \
	-out_scaf data/mlsmr_hscaf_scaf_dev.smi \
	-stereo \
	-inc_mol -inc_scaf \
	-dbdir scratch
#
#
