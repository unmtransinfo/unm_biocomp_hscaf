#!/bin/sh
#
set -x
#
#
if [ `uname -s` = "Darwin" ]; then
	appdir=/Users/app
else
	appdir=/home/app
fi
#
export JAVA_HOME=/usr/lib64/jvm/java
#export JAVA_HOME=/usr/lib64/jvm/java-1.5.0
export JAVA=$JAVA_HOME/bin/java
#
#INFILE=$HOME/chemaxon/hscaf/data/testset.smi
#INFILE=$HOME/Download/mlsmr_500.smi
#INFILE=$HOME/Download/mlsmr_700.smi
#INFILE=$HOME/Download/mlsmr_1000.smi
INFILE=$HOME/Download/mlsmr_5000.smi
#
# (Start $JAVA_HOME/bin/jvisualvm)
#
CLASSPATH=$HOME/src/java/lib/unm_biocomp_hscaf.jar
#CLASSPATH=$appdir/lib/unm_biocomp_hscaf.jar
CLASSPATH=$CLASSPATH:$HOME/src/java/lib/unm_biocomp_util.jar
CLASSPATH=$CLASSPATH:$appdir/lib/jchem.jar
CLASSPATH=$CLASSPATH:$appdir/lib/berkeleydb.jar
#
$JAVA \
	-classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds \
	-i $INFILE \
	-out_scaf $HOME/chemaxon/hscaf/data/z_scaf.smi \
	-v
#
