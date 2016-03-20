#!/bin/bash

CMDFILE=/usr/local/bin/datalog.cmd
WORKDIR=/tmp/stable-stats
SCPDEST="jra@meow:/var/www/febo.com/pages/plots/chronos/"
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
	# run stable-stats program
	stable-stats.pl -t $TAU -n $ALIGN -g $TITLE -s $SUB -b $BASE -f $DFILE

done < $CMDFILE
IFS=$OLDIFS

# move files to web server
scp $WORKDIR/*.html $SCPDEST
scp $WORKDIR/*.png $SCPDEST/images/
rm $WORKDIR/*
