#!/bin/sh
#
set -x
#
LIBDIR=$HOME/src/java/lib
CLASSPATH=$LIBDIR/unm_biocomp_hscaf.jar
CLASSPATH=$CLASSPATH:$LIBDIR/unm_biocomp_db.jar
CLASSPATH=$CLASSPATH:$LIBDIR/unm_biocomp_util.jar
LIBDIR=/home/app/lib
CLASSPATH=$CLASSPATH:$LIBDIR/jchem.jar
CLASSPATH=$CLASSPATH:$LIBDIR/berkeleydb.jar
#
IFILE=$HOME/projects/badapple/data/pubchem_mlsmr.smi
NMAX=50000
#
thisdir=`dirname $0`
#
JAVA_OPTS='-Xmx4G -Xms4G'
#
java $JAVA_OPTS -classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds \
	-i $IFILE \
	-o $thisdir/data/z_hscaf_out_rdb.smi \
	-out_scaf $thisdir/data/z_hscaf_scaf_rdb.smi \
	-inc_mol \
	-inc_scaf \
	-inc_link \
	-inc_chain \
	-nmax $NMAX \
	-maxmol 80 \
	-maxrings 5 \
	-scaflist_append2title \
	-rdb \
	-rdb_host "localhost" \
	-rdb_port "5432" \
	-rdb_name "hscafscratch" \
	-rdb_schema "public" \
	-rdb_user "jjyang" \
	-rdb_pw "assword" \
	-rdb_predelete \
	-rdb_keep \
	-v
#
#	-rdb_predelete \
#	-rdb_keep \
#
