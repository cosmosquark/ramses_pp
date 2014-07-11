#!/bin/bash
SCRIPT="`readlink -e $0`"
SCRIPTPATH="`dirname $SCRIPT`"
DIRPATH="$(dirname $SCRIPTPATH)"
cp config.py.in config.py
sed -i 's,'"..ramses_f90"','"${DIRPATH}/ramses_f90"',g' config.py
sed -i 's,'"..data"','"${DIRPATH}/data"',g' config.py
sed -i 's,'"..simulations"','"${DIRPATH}/simulations"',g' config.py
sed -i 's,'"..applications"','"${DIRPATH}/applications"',g' config.py
sed -i 's,'"..topdog"','"${DIRPATH}"',g' config.py
sed -i 's,'"..thing"','"${DIRPATH}/ramses_pp"',g' config.py
sed -i 's,'"..really"','"${DIRPATH}/ramses_pp/applications"',g' config.py
