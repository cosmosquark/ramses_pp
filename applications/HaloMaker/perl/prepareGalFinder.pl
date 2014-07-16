#!/usr/bin/perl
use strict;
use warnings;

##
## THIS WOULD WORK IF GalFinder DID ... BUT IT DOES NOT ... ## 
##

# define list of snapshots to process
my @snapList = (9,17);
#for (my $k = 10-9; $k<128-9; $k++) {$snapList[$k] = $k+10;}

# define snapshot directory (dir containing all output_xxxxx directories)
my $snapDir = "/data/blaizot/TMP/Paulina/";

# define (and create if needed) global output directory (which will host many sub-directories...)
my $runDir = "/data/blaizot/TMP/Paulina/Catalogs/Gal/";
if (!-e $runDir) {system("mkdir $runDir");}

# copy code over to run directory, as a documentation ... 
if (!-e $runDir."/code") {system("mkdir $runDir"."/code");}
system("cp -r ../\* $runDir"."/code/.");

# define executable file as the above copied executable (so we can link with no risk). 
my $execFile  = "$runDir"."/code/f90/GalFinder";
my $cmdFile   = "$runDir" . "cmdFile.sh";    # file which will contain instructions to run HaloFinder


# LOOP ON SNAPSHOTS 
open(CFILE,">$cmdFile");
foreach (@snapList)
  {
    # convert snapshot number into i5.5 format (e.g. 00001)
    my $i   = $_;
    my $num = "$i";
    if ($i < 10000) {$num = "0"."$num";}
    if ($i < 1000)  {$num = "0"."$num";}
    if ($i < 100)   {$num = "0"."$num";}
    if ($i < 10)    {$num = "0"."$num";}

    # define location of simulation files (snapshots)
    my $fileRoot  = "$snapDir" . "/output_$num/";

    # create working directory for snapshot
    my $dir = "$runDir"."$i";
    my $cmd = "mkdir $dir";
    system($cmd);
    
    # write input_HaloMaker.dat file
    write_input_hm($dir);

    # link executable and stuff to running dir.
    $cmd = "ln -s $execFile $dir"; 
    system($cmd);

    # write inputfiles_HaloMaker.dat file
    my $file = "$dir"."/inputfiles_HaloMaker.dat";
    open(JFILE,">$file");
    print JFILE "\'$fileRoot \'  Ra3   1  $num \n";
    close(JFILE);
    
    print CFILE "cd $dir && ./GalFinder > log.out \n";
  }
close(CFILE);


sub write_input_hm {
    my $filename  = "$_[0]"."/input_HaloMaker.dat";
    open(FILE,">$filename");
    # MN PARAMETERS 
    print FILE "af              = 1.0        	! Today expansion factor == 1+zi \n";
    print FILE "lbox            = 28.5714    	! Size of the box at z=0 (NB: in Mpc! ==> 479 Mpc/h=684.286 Mpc) \n";
    print FILE "H_f             = 70.000      	! Hubble constant (in Km/s/Mpc) at z=0 \n";
    print FILE "omega_f         = 0.3        	! Omega matter at z=0 (ex: 0.3) \n";
    print FILE "lambda_f        = 0.7        	! Omega lambda at z=0 (ex: 0.7) \n";    
    print FILE "FlagPeriod      = 1           	! pericidicity flag choose 1 for periodic boundary conditions 0 else \n";
    print FILE "npart           = 20            ! Minimum number of particles per halo (ex:10) \n";
    print FILE "method          = MSM        	! choose between FOF, HOP, DPM, MHM (view ReadMe for details). \n";
    print FILE "cdm             = .false.    	! use center of mass instead of densest particle. \n";
    print FILE "b               = 0.2        	! fof parameter b (usually b=0.2) \n";
    print FILE "nsteps          = 1         	! Number of time steps to analyse \n";
    # adaptaHOP parameters
    print FILE "nvoisins        = 20        	! parameter for adaptahop (usually 20)\n";
    print FILE "nhop            = 20         	! parameter for adaptahop (usually 20)\n";
    print FILE "rhot            = 2000.       	! parameter for adaptahop (80 coreespond to b = 0.2)\n";
    print FILE "fudge           = 5.          	! parameter for adaptahop (usually 4.) \n";
    print FILE "fudgepsilon     = 1.d-6         ! parameter for adaptahop (usually 0.05)\n";
    print FILE "alphap          = 1.          	! parameter for adaptahop (usually 1.)\n";
    print FILE "verbose         = .true.     	! verbose parameter for both halo finder\n";
    print FILE "megaverbose     = .false.     	! parameter for adaptahop\n";
    close(FILE);
}


