#!/bin/sh
#
set -x
#
JAVA=java
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
## By parts...
#
INFILE=data/mlp_1000.smi
PREFIX=data/mlp_1000
SCRATCHDIR=scratch
#
if [ `ls $SCRATCHDIR` ]; then
	echo "ERROR: scratch dir \"$SCRATCHDIR\" not empty."
	exit
fi
#
NMOL=`cat $INFILE |wc -l`
NJOB=100
#
nskip=0
ndone=0
i=0
ntodo=$NMOL
#
while [ "$ntodo" -gt 0 -a "$nskip" -lt $NMOL ]; do
	i=`expr $i + 1`
	i=`printf "%03d" $i`
	#
	#opts="-v"
	opts=""
	if [ $i -gt 1 ]; then
		nskip=`expr $nskip + $NJOB`
		opts=" -resume -nskip $nskip"
	fi
	if [ $NJOB -le $ntodo ]; then
		opts="$opts -keepdb -nmax $NJOB"
	fi
	#
	$JAVA -classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds \
		$opts \
		-i $INFILE \
		-o ${PREFIX}_byparts_out_${i}.smi \
		-out_scaf ${PREFIX}_byparts_hscaf_scaf_${i}.smi \
		-inc_mol \
		-inc_scaf \
		-inc_link \
		-inc_chain \
		-scaflist_append2title \
		-dbdir $SCRATCHDIR
	#
	ndone_this=`cat ${PREFIX}_byparts_out_${i}.smi |wc -l`
	ndone=`expr $ndone + $ndone_this`
	ntodo=`expr $NMOL - $ndone`
done
#
#
