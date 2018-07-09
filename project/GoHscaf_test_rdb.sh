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
#
thisdir=`dirname $0`
#
java -classpath $CLASSPATH edu.unm.health.biocomp.hscaf.hier_scaffolds \
	-i $HOME/data/smi/hscaf_testset.smi \
	-o $thisdir/data/z_hscaf_out_rdb.smi \
	-out_scaf $thisdir/data/z_hscaf_scaf_rdb.smi \
	-inc_mol \
	-inc_scaf \
	-inc_link \
	-inc_chain \
	-scaflist_append2title \
	-rdb \
	-rdb_host "localhost" \
	-rdb_port "5432" \
	-rdb_name "hscafscratch" \
	-rdb_schema "public" \
	-rdb_user "jjyang" \
	-rdb_pw "assword" \
	-rdb_predelete \
	-vvv
#
#	-rdb_predelete \
#	-rdb_keep \
#
