#!/bin/sh
#
set -x
#
CLASSPATH=/home/app/lib/unm_biocomp_hscaf.jar
CLASSPATH=$CLASSPATH:/home/app/ChemAxon/JChem/lib/jchem.jar
CLASSPATH=$CLASSPATH:/usr/share/java/berkeleydb.jar
#
java -classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds_common
#
#	-i $HOME/Download/mlsmr.smi \
#	-o data/mlsmr_hscaf_out.smi \
#	-out_scaf data/mlsmr_hscaf_scaf.smi \
#	-scaflist_append2title \
#	-scratchdir scratch \
#	-v
#
#
