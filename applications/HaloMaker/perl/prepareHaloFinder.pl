#!/Usr/bin/perl
use strict;
use warnings;

# define list of snapshots to process
my @snapList = (1);
for (my $k = 0; $k<238; $k++) {$snapList[$k] = $k+1;}

# define snapshot directory (dir containing all output_xxxxx directories)
#my $snapDir = "/data/blaizot/TMP/Paulina/";
my $snapDir = "/sps/cral/leo/Ramses_v3.07/boxlen20_n1024_W5/Zoom331/Feedbk/";

# define (and create if needed) global output directory (which will host many sub-directories...)
#my $runDir = "/data/blaizot/TMP/Paulina/Catalogs/DM/";
my $runDir = "/sps/cral/leo/Ramses_v3.07/boxlen20_n1024_W5/Zoom331/Feedbk/Catalogs/HaloMaker/";
if (!-e $runDir) {system("mkdir $runDir");}

# copy code over to run directory, as a documentation ... 
if (!-e $runDir."/code") {system("mkdir $runDir"."/code");}
system("cp -r ../f90/\* $runDir"."/code/.");

# define executable file as the above copied executable (so we can link with no risk). 
my $execFile  = "$runDir"."/code/HaloFinder";
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
    #write_input_hm($dir);
    write_input_hm_sub($dir);

    # link executable and stuff to running dir.
    $cmd = "ln -s $execFile $dir"; 
    system($cmd);

    # write inputfiles_HaloMaker.dat file
    my $file = "$dir"."/inputfiles_HaloMaker.dat";
    open(JFILE,">$file");
    print JFILE "\'$fileRoot \'  Ra3   1  $num \n";
    close(JFILE);
    
    print CFILE "cd $dir && ./HaloFinder > log.out \n";
  }
close(CFILE);

## write job file -> plateform dependant -> ccage openmpi
#open (JFILE,">rats.job");
#print JFILE "#!/usr/local/bin/bash -l \n";
#print JFILE "#\$ -P P_cral \n";
#print JFILE "#\$ -l sps=1 \n";
#print JFILE "#\$ -l ct=18000       # cpu in s \n";
#print JFILE "#\$ -l vmem=3G        # (required memory) \n";
#print JFILE "#\$ -l fsize=2G       # required DISK \n";
#print JFILE "#\$ -N rats \n";
#print JFILE "#\$ -o log.txt \n";
#print JFILE "#\$ -e errors.txt \n";
#print JFILE "#\$ -cwd              # run on current working directory \n";
#print JFILE "#\$ -pe openmpi 3 \n";
#print JFILE "#\$ -q pa_medium \n";
#print JFILE "echo '   job started at' \n";
#print JFILE "date \n";
#print JFILE "cd $runDir \n";
#print JFILE ". /usr/local/shared/bin/openmpi_env.sh \n";
#print JFILE "source /usr/local/intel/ifort/bin/ifortvars.sh \n";
#print JFILE "/usr/local/openmpi/bin/mpiexec --mca btl ^udapl,openib --mca btl_tcp_if_include eth2 -n \$NSLOTS /sps/cral/leo/RATS_vL#eo/MultiProcStarter/mpp_starter $cmdFile > run.log \n";
#print JFILE "echo '   job ended at' \n";
#print JFILE "date \n";
#close(JFILE);

# write job file -> plateform dependant -> ccage mpich2 (doesn't need high speed communication)
my $file = "$runDir"."/rats.job";
open (JFILE,">$file");
print JFILE "#!/usr/local/bin/bash -l \n";
print JFILE "#\$ -P P_cral \n";
print JFILE "#\$ -l sps=1 \n";
print JFILE "#\$ -l ct=44000       # cpu in s \n";
print JFILE "#\$ -l vmem=4G        # (required memory) \n";
print JFILE "#\$ -l fsize=2G       # required DISK \n";
print JFILE "#\$ -N rats \n";
print JFILE "#\$ -o log.txt \n";
print JFILE "#\$ -e errors.txt \n";
print JFILE "#\$ -cwd              # run on current working directory \n";
print JFILE "#\$ -pe mpich2_hf 9 \n";
#print JFILE "#\$ -pe mpich2 3 \n";
print JFILE "#\$ -q pa_long \n";
print JFILE "echo '   job started at' \n";
print JFILE "date \n";
print JFILE "cd $runDir \n";
print JFILE ". /usr/local/shared/bin/mpich2_hydra_env.sh \n";
print JFILE "export MPICH2_ROOT=/usr/local/mpich2/hydra \n";
print JFILE "/usr/local/mpich2/hydra/bin/mpiexec -rmk sge -iface eth2 -np \$NSLOTS /sps/cral/leo/RATS_vLeo/MultiProcStarter/mpp_starter $cmdFile > run.log \n";
#print JFILE ". /usr/local/shared/bin/mpich2_env.sh \n";
#print JFILE "export MPD_CON_EXT='sge_\$JOB_ID.\$SGE_TASK_ID'  \n";
#print JFILE "export MPICH2_ROOT=/usr/local/mpich2/ \n";
#print JFILE "/usr/local/mpich2/bin/mpiexec -machinefile \$TMPDIR/machines -n \$NSLOTS /sps/cral/leo/RATS_vLeo/MultiProcStarter/mpp_starter $cmdFile > run.log \n";
print JFILE "echo '   job ended at' \n";
print JFILE "date \n";
close(JFILE);


#open(JFILE,">toto.job");
#print JFILE "#!/usr/local/bin/bash \n";
#print JFILE "#PBS -l platform=LINUX \n";
#print JFILE "#PBS -l proc=6 \n";
#print JFILE "#PBS -l u_sps_cral=6 \n";
#print JFILE "#PBS -l machine=2 \n";
#print JFILE "#PBS -l ptype=MPICH2 \n";
#print JFILE "#PBS -l M=4GB \n";
#print JFILE "#PBS -l T=4000000 \n";
#print JFILE "#PBS -N HM \n";
#print JFILE "#PBS -o log.txt \n";
#print JFILE "#PBS -e errors.txt \n";
#print JFILE " \n";
#print JFILE "cd \$PBS_O_WORKDIR \n";
#print JFILE ". /usr/local/shared/bin/mpich2_env.sh \n";
#print JFILE "mpdboot -n \$BQS_HOSTNUMBER -f \$BQS_HOSTLISTPATH.GIGA --ifhn=\$(hostname)p --verbose --ncpus=\${BQS_PROCBYHOST} --rsh=/usr/local/products/bqs/bqsrsh \n"#;
#print JFILE "mpiexec -machinefile \$BQS_PROCLISTPATH.GIGA.MPIEXEC -np \$BQS_PROCNUMBER ../../MultiProcStarter/mpp_starter $cmdFile > run.log \n";
#print JFILE "mpdallexit \n";
#print JFILE " \n";
#close(JFILE);



sub write_input_hm {
    my $filename  = "$_[0]"."/input_HaloMaker.dat";
    open(FILE,">$filename");
    # MN PARAMETERS 
    print FILE "af              = 1.0        	! Today expansion factor == 1+zi \n";
    print FILE "lbox            = 28.5714    	! Size of the box at z=0 (NB: in Mpc! ==> 479 Mpc/h=684.286 Mpc) \n";
    print FILE "H_f             = 70.000      	! Hubble constant (in Km/s/Mpc) at z=0 \n";
    print FILE "omega_f         = 0.277        	! Omega matter at z=0 (ex: 0.3) \n";
    print FILE "lambda_f        = 0.723        	! Omega lambda at z=0 (ex: 0.7) \n";    
    print FILE "FlagPeriod      = 1           	! pericidicity flag choose 1 for periodic boundary conditions 0 else \n";
    print FILE "npart           = 20            ! Minimum number of particles per halo (ex:10) \n";
    print FILE "method          = MSM        	! choose between FOF, HOP, DPM, MHM (view ReadMe for details). \n";
    print FILE "cdm             = .false.    	! use center of mass instead of densest particle. \n";
    print FILE "b               = 0.2        	! fof parameter b (usually b=0.2) \n";
    print FILE "nsteps          = 1         	! Number of time steps to analyse \n";
    # adaptaHOP parameters
    print FILE "nvoisins        = 20        	! parameter for adaptahop (usually 20)\n";
    print FILE "nhop            = 20         	! parameter for adaptahop (usually 20)\n";
    print FILE "rhot            = 80.       	! parameter for adaptahop (80 coreespond to b = 0.2)\n";
    print FILE "fudge           = 5.          	! parameter for adaptahop (usually 4.) \n";
    print FILE "fudgepsilon     = 1.d-6         ! parameter for adaptahop (usually 0.05)\n";
    print FILE "alphap          = 1.          	! parameter for adaptahop (usually 1.)\n";
    print FILE "verbose         = .true.     	! verbose parameter for both halo finder\n";
    print FILE "megaverbose     = .false.     	! parameter for adaptahop\n";
    close(FILE);
}

sub write_input_hm_sub {
    my $filename  = "$_[0]"."/input_HaloMaker.dat";
    open(FILE,">$filename");
    # MN PARAMETERS 
    print FILE "af              = 1.0        	! Today expansion factor == 1+zi \n";
    print FILE "lbox            = 28.5714    	! Size of the box at z=0 (NB: in Mpc! ==> 479 Mpc/h=684.286 Mpc) \n";
    print FILE "H_f             = 70.000      	! Hubble constant (in Km/s/Mpc) at z=0 \n";
    print FILE "omega_f         = 0.277        	! Omega matter at z=0 (ex: 0.3) \n";
    print FILE "lambda_f        = 0.723        	! Omega lambda at z=0 (ex: 0.7) \n";    
    print FILE "FlagPeriod      = 1           	! pericidicity flag choose 1 for periodic boundary conditions 0 else \n";
    print FILE "npart           = 20            ! Minimum number of particles per halo (ex:10) \n";
    print FILE "method          = MSM        	! choose between FOF, HOP, DPM, MHM (view ReadMe for details). \n";
    print FILE "cdm             = .false.    	! use center of mass instead of densest particle. \n";
    print FILE "b               = 0.2        	! fof parameter b (usually b=0.2) \n";
    print FILE "nsteps          = 1         	! Number of time steps to analyse \n";
    # adaptaHOP parameters
    print FILE "nvoisins        = 20        	! parameter for adaptahop (usually 20)\n";
    print FILE "nhop            = 20         	! parameter for adaptahop (usually 20)\n";
    print FILE "rhot            = 80.       	! parameter for adaptahop (80 coreespond to b = 0.2)\n";
    print FILE "fudge           = 5.          	! parameter for adaptahop (usually 4.) \n";
    print FILE "fudgepsilon     = 1.d-6         ! parameter for adaptahop (usually 0.05)\n";
    print FILE "alphap          = 1.          	! parameter for adaptahop (usually 1.)\n";
    print FILE "verbose         = .true.     	! verbose parameter for both halo finder\n";
    print FILE "megaverbose     = .false.     	! parameter for adaptahop\n";
    # subbox parameters
    print FILE "xmin            = -0.1 \n";
    print FILE "xmax            = 0.1 \n";
    print FILE "ymin            = -0.1 \n";
    print FILE "ymax            = 0.1 \n";
    print FILE "zmin            = -0.1 \n";
    print FILE "zmax            = 0.1 \n";
    print FILE "padlength       = 0.05 \n";
    print FILE "DangerPadLength = 0. \n";

    close(FILE);
}


