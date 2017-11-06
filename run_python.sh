#!/bin/sh 
# just a simple script to run python in a localised enviroment... i.e here
# written by B B Thompson (@Cosmosquark)


## this program is used to run python in the localised enviroment set by ramses_pp. This will save you from configuring python paths and allows you to run ramses_pp as a seperate entity
##############################

# finds the current location of this script, please do not change this location
SCRIPT="`readlink -e $0`"
SCRIPTPATH="`dirname $SCRIPT`"
PYTHONLOC="/python"

# runs ramses_pp from the directory above (e.g so we can use from ramses_pp import config for ramses_pp specific things)

echo "Location of ramses_pp"
echo $SCRIPTPATH
PYTHONPATH=$SCRIPTPATH$PYTHONLOC

echo $PYTHONPATH
# now we want parent directory of this script, i.e the ramses_pp root name

#parentdir="$(dirname $SCRIPTPATH)"

#echo $parentdir

#PYTHONPATH=$parentdir



export PYTHONPATH

/usr/bin/python $@
