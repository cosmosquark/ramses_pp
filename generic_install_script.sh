#!/bin/bash

#installdir="numpy"
#installdir="scipy"
#installdir="Cython"
#installdir="pynbody"
installdir="pymses"
fcompiler="gfortran"


SCRIPT="`readlink -e $0`"
SCRIPTPATH="`dirname $SCRIPT`"
PYTHONLOC="/python"

# runs ramses_pp from the directory above (e.g so we can use from ramses_pp import config for ramses_pp specific things)

echo "Location of ramses_pp"
echo $SCRIPTPATH
PYTHONPATH=$SCRIPTPATH$PYTHONLOC
export PYTHONPATH

cd python
cd $installdir

INSTALLOC="`pwd`"


# numpy installation command
#/usr/bin/python $@ setup.py build --fcompiler=$fcompiler

# scipy installation command
/usr/bin/python $@ setup.py install --prefix=$INSTALLOC
