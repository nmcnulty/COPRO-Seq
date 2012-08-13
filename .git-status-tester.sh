#!/bin/sh
# Author: Nick Semenkovich <semenko@alum.mit.edu>
# License: MIT
# Version: 1.1 [2012.08.13]

# Move to our directory, if from a relative path.
# TODO: This readlink breaks on BSD/OSX.
SCRIPT=$(readlink -f "$0")
SCRIPTPATH=`dirname "$SCRIPT"`
cd $SCRIPTPATH

git remote update >/dev/null 2>&1

gitstatus=`git log master..origin/master`

if [ "$gitstatus" ]
then
    if tty -s
    then
	USER_TTY=`tty`
	echo "\n\a\t\t@@@ WARNING @@@ - This software is outdated." > $USER_TTY
	echo "\n\nThere is an updated version in git with these changes:\n" > $USER_TTY
	echo "$gitstatus" > $USER_TTY
	echo "\n\n\t\t@@@ Please run: git pull @@@\n" > $USER_TTY
	sleep 2
    fi
fi

#END