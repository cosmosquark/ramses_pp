#!/bin/sh 
# just a simple script to run python in a localised enviroment... i.e here
# written by B B Thompson (@Cosmosquark)


SCRIPT="`readlink -e $0`"
SCRIPTPATH="`dirname $SCRIPT`"

echo "Location of ramses_pp"
echo $SCRIPTPATH

PYTHONPATH=$SCRIPTPATH
export PYTHONPATH

/usr/bin/python $@
