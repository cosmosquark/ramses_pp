#!/usr/bin/perl
use strict;
use warnings;

# define list of snapshots to process
my @snapList = (10);
for (my $k = 10-9; $k<128-9; $k++) {$snapList[$k] = $k+10;}

# define (and create if needed) global output directory (which will host many sub-directories...)
#my $runDir    = "/data/blaizot/Horizon/1024-100h-1Mpc-W3/HaloMaker/MSM/5x5x5/";
#my $runDir    = "/data/blaizot/TEST_HM/512-100h-1Mpc-W3/TEST/4x4x4/";
my $runDir = "/data33/blaizot/MN/HaloMaker/3x3x3/";
#my $runDir    = "/data/blaizot/TEST_HM/MN/MSM/test/";
if (!-e $runDir) {system("mkdir $runDir");}

# copy code over to run directory, as a documentation ... 
if (!-e $runDir."/code") {system("mkdir $runDir"."/code");}
system("cp \* $runDir"."/code/.");

# define executable file as the above copied executable (so we can link with no risk). 
my $execFile  = "$runDir"."/code/HaloMaker";
my $cmdFile   = "$runDir" . "cmdFile.sh";    # instructions to run HM on all sub-boxes
my $mergeFile = "$runDir" . "mergeFile.sh";  # instructions to merge sub-box files

# define sub-boxes properties 
my $nsbx = 3;
my $sbdx = 1.0 / $nsbx;
my $nsby = 3;
my $sbdy = 1.0 / $nsby;
my $nsbz = 3;
my $sbdz = 1.0 / $nsbz;

# write out details of run into file used to merge catalogues 
open(MFILE,">$mergeFile");
my $nsb = $nsbx * $nsby * $nsbz;
print MFILE scalar @snapList . " $nsb \n";
foreach (@snapList)
  {
    print MFILE "$_ \n";
  }
close(MFILE);

open(CFILE,">$cmdFile");
# LOOP ON SNAPSHOTS
foreach (@snapList)
  {
    my $i   = $_;
    my $num = "$i";
    if ($i < 10000) {$num = "0"."$num";}
    if ($i < 1000)  {$num = "0"."$num";}
    if ($i < 100)   {$num = "0"."$num";}
    if ($i < 10)    {$num = "0"."$num";}

    # MareNostrum stuff : get FudgeEpsilon from a file ... 
    my $fudgeEpsilon = get_fudge($i);

    # define location of simulation files (snapshots)
    #my $fileRoot  = "/simulations/Horizon/1024-100h-1Mpc-W3/Snapshots/snapshot_$num"."/snapshot_";
    #my $fileRoot  = "/simulations/Horizon/512-100h-1Mpc-W3/Snapshots/snapshot_";
    #my $fileRoot  = "/simulations/MN/outputs/output_00038/";
    my $fileRoot  = "/data33/blaizot/MN_SNAPS/output_$num/";

    my $dir = "$runDir"."$i";
    my $cmd = "mkdir $dir";
    system($cmd);
    
    # LOOP ON SUBBOXES
    my $subboxnum = 0;
    for (my $ix = 0; $ix < $nsbx; $ix++) 
    {
	for (my $iy = 0; $iy < $nsby; $iy++) 
	{
	    for (my $iz = 0; $iz < $nsbz; $iz++) 
	    {
		# create directory for subbox run
		my $sbdir = "$dir" . "/$subboxnum";
		$cmd = "mkdir $sbdir";
		system($cmd);
		# define sub-box boundaries
		my $xmin = $ix * $sbdx - 0.5;
		my $xmax = ($ix + 1) * $sbdx - 0.5;
		my $ymin = $iy * $sbdy - 0.5;
		my $ymax = ($iy + 1) * $sbdy - 0.5;
		my $zmin = $iz * $sbdz - 0.5;
		my $zmax = ($iz + 1) * $sbdz - 0.5;
		# write input_HaloMaker.dat file
		write_input_hm($sbdir,$subboxnum,$xmin,$xmax,$ymin,$ymax,$zmin,$zmax,$fudgeEpsilon);
		# copy executable and stuff to running dir.
		$cmd = "ln -s $execFile $sbdir"; 
		system($cmd);
		# write inputfiles_HaloMaker.dat file
		my $file = "$sbdir"."/inputfiles_HaloMaker.dat";
		open(JFILE,">$file");
		#print JFILE "\'$fileRoot"."$num \'  Gd   1  $num \n";
		print JFILE "\'$fileRoot \'  Ra   1  $num \n";
		close(JFILE);
		print CFILE "cd $sbdir && ./HaloMaker > log.out \n";
		$subboxnum++ ;
	    }
	}
    }
}
close(CFILE);


sub write_input_hm {
    my $filename  = "$_[0]"."/input_HaloMaker.dat";
    my $subboxnum = $_[1];
    my $xmin      = $_[2];
    my $xmax      = $_[3];
    my $ymin      = $_[4];
    my $ymax      = $_[5];
    my $zmin      = $_[6];
    my $zmax      = $_[7];
    my $fudgeEps  = $_[8];
    open(FILE,">$filename");
    # MN PARAMETERS 
    print FILE "af             =  1.0        	! Today expansion factor == 1+zi \n";
    print FILE "lbox           =  71.428    	! Size of the box at z=0 (NB: in Mpc! ==> 479 Mpc/h=684.286 Mpc) \n";
    print FILE "H_f            =  70.000      	! Hubble constant (in Km/s/Mpc) at z=0 \n";
    print FILE "omega_f        =  0.3        	! Omega matter at z=0 (ex: 0.3) \n";
    print FILE "lambda_f       =  0.7        	! Omega lambda at z=0 (ex: 0.7) \n";    
    # 256-100h-1Mpc-W3 parameters
    #print FILE "af              =  1.0        	! Today expansion factor == 1+zi \n";
    #print FILE "lbox            =  136.9863    	! Size of the box at z=0 (NB: in Mpc! ==> 479 Mpc/h=684.286 Mpc) \n";
    #print FILE "H_f             =  73.000      	! Hubble constant (in Km/s/Mpc) at z=0 \n";
    #print FILE "omega_f         =  0.24        	! Omega matter at z=0 (ex: 0.3) \n";
    #print FILE "lambda_f        =  0.76        	! Omega lambda at z=0 (ex: 0.7) \n";
    print FILE "FlagPeriod      =  1           	! pericidicity flag choose 1 for periodic boundary conditions 0 else \n";
    print FILE "npart           =  20            ! Minimum number of particles per halo (ex:10) \n";
    print FILE "method          =  MSM        	! choose between FOF, HOP, DPM, MHM (view ReadMe for details). \n";
    print FILE "cdm             =  .false.    	! use center of mass instead of densest particle. \n";
    print FILE "b               =  0.2        	! fof parameter b (usually b=0.2) \n";
    print FILE "nsteps          =  1         	! Number of time steps to analyse \n";
    # jeje parameters
    print FILE "subboxnum       = $subboxnum \n";
    print FILE "xmin            = $xmin          ! define sub-box limits (in code units) \n";
    print FILE "xmax            = $xmax          ! -> the search will extend from xmin-padlehgth (below) \n";
    print FILE "ymin            = $ymin          !    to xmax+padlength, but only halos within ]xmin,xmax] \n";
    print FILE "ymax            = $ymax          !    are kept. \n";
    print FILE "zmin            = $zmin \n";
    print FILE "zmax            = $zmax \n";
    # for Lbox = 100Mpc/h:
    #print FILE "padlength       = 0.05           ! buffer size (in comoving box units) \n";
    #print FILE "DangerPadLength = 0.01         ! danger zone within the par (same units as padlength) \n";
    # for Lbox = 50Mpc/h: (try and see ... )
    print FILE "padlength       = 0.08          ! buffer size (in comoving box units) \n";
    print FILE "DangerPadLength = 0.007         ! danger zone within the par (same units as padlength) \n";
    # adaptaHOP parameters
    print FILE "nvoisins        =  20        	! parameter for adaptahop (usually 20)\n";
    print FILE "nhop            =  20         	! parameter for adaptahop (usually 20)\n";
    print FILE "rhot            =  80.       	! parameter for adaptahop (80 coreespond to b = 0.2)\n";
    print FILE "fudge           = 5.          	! parameter for adaptahop (usually 4.) \n";
    print FILE "fudgepsilon     = $fudgeEps   	! parameter for adaptahop (usually 0.05)\n";
    print FILE "alphap          = 1.          	! parameter for adaptahop (usually 1.)\n";
    print FILE "verbose         = .true.     	! verbose parameter for both halo finder\n";
    print FILE "megaverbose     = .false.     	! parameter for adaptahop\n";
    close(FILE);
}

sub get_fudge{
    my $wantedsnap=$_[0];
    my $fudgeFile = "/data33/blaizot/MN_SNAPS/snapshot_fudges.txt";
    open FFILE, "<$fudgeFile";
    while (my $line = <FFILE>){
	chomp($line);
	my $isnap;
	my $fue;
	($isnap,$fue) = split(/,/,$line);
	$isnap =~ s/^\s+//; # remove leading spaces
	$isnap =~ s/\s+$//; # remove trailing spaces
	$fue   =~ s/^\s+//; # remove leading spaces
	$fue   =~ s/\s+$//; # remove trailing spaces
	if ($isnap == $wantedsnap) {return($fue);}
    }
    close(FFILE);
}

