

Installing PNG needs impi loaded.. same with freetype and the earlier stuff.
SQLite wants impi unloaded... wants the other intel MPI libraries
numpy wants impi loaded. ... ordering matters... impi/openmpi needs to be loaded AFTER openmpi/1.6.5/gcc (i.e a higher number).. but openmpi/1.6.5/gcc is needed for the installation of distribute
h5py wants impi unloaded... wants the other intel libraries... which seems to be working for everything else
be sure to load impi later

modules needd

  1) icomp/11.1.073             3) imkl/10.2.5.035            5) openmpi/1.6.5/gcc          7) scalapack/2.0.2
  2) openmpi/1.4.4-intel_11.1   4) lapack/3.5.0               6) python/2.7.6               8) ffmpeg/0.6.1
  9) impi/4.0.3

also, you may need to do this

unset CFLAGS

other things needed

0MQ does not easily install.. but is not needed.. you can later address this via
pip install pyzmq

check installation

import numpy
import scipy
import matplotlib
import cython
import h5py

  openmpi/1.4.4-intel_11.1 may now be the best one to be "last" out of the series of modules

pynbody installation

try easy_install pynbody or pip install pynbody
failing that.. try and install it manually outside the enviroment.. and move the build into the ramses_pp directory.. or python site package directory

pymses installation (requirements wxpython, numexpr, tables)

easy_install wxpython
(note.. wxpython instal scripts seem to be missing... if this is the case.. then
wget http://downloads.sourceforge.net/wxpython/wxPython-src-3.0.0.0.tar.bz2
tar -jxvf wxPython-src-3.0.0.0.tar.bz2
cd wxPython-src-3.0.0.0

and try and install that...

failing that

easy_install pil
will do for most functionality
)

easy_install numexpr
cd /yt-x86_64/src/hdf5-1.8.11
(if the hdf5 dir does not exist.. run ./configure then make)
cd hdf5
export HDF5_DIR=`pwd`
easy_install tables 

: Parent module 'numpy.distutils' not found while handling absolute import
  from numpy.distutils import log

ignore this error


finally run

ramses_pp/gen_config.sh

and we're ready to go
