#!/bin/bash
SCRIPT="`readlink -e $0`"
SCRIPTPATH="`dirname $SCRIPT`"
cp config.py.in config.py
sed -i 's,'"..ramses_f90"','"${SCRIPTPATH}/ramses_f90"',g' config.py
sed -i 's,'"..data"','"${SCRIPTPATH}/data"',g' config.py
