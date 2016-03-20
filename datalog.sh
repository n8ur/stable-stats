#!/bin/bash

CMDFILE=/usr/local/bin/datalog.cmd
WORKDIR=/tmp/stable-stats
SCPDEST="jra@meow:/var/www/febo.com/pages/plots/chronos/"
THERMDAT=/data/therm.dat
TAU=3600
ALIGN="center"

mkdir -p $WORKDIR
cd $WORKDIR

OLDIFS=$IFS
IFS="|"
[ ! -f $CMDFILE ] &while read DFILE DUT_CH REF_CH TITLE SUB COMMENT
do

	# skip blank lines and lines starting with #
	case $DFILE in ''|\#*) continue ;; esac
	# get base
	BASE=$(basename "$DFILE" | cut -d. -f1)
	# create full comment string
	FULL_COMMENT="$TITLE $COMMENT"

	# run datalogger program
	/usr/local/bin/hp5334-logger.pl -a $DUT_CH -b $REF_CH -s 10 \
		-d $THERMDAT -c $FULL_COMMENT -z ab -w -f $DFILE


done < $CMDFILE
IFS=$OLDIFS
