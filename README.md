Collection of post-processing scripts for RAMSES  using Pymses, Pynbody and YT

Dependencies (as of 10/04/14): YT, Pymses, Pynbody, mpi4py, numpy, matplotlib
Optional (recommended) dependencies (as of 10/04/14): MPlayer (mencoder), ffmpeg

Install yt3.0 via the installation script. As of 10/04/14, HOP will not work with RAMSES output. To fix this,
cd to '/$PATHTOYT/yt-X/src/yt-hg/' and run 'hg up development' followed by 'python setup.py develop'.

Once completed, sym-link ramses_pp to the directory '/$PATHTOYT/yt-X/lib/python-2.7/site-packages/'.

Job done!

Developed by Ben Thompson (Cosmosquark) founded by David Sullivan (daibandofullpelt)
