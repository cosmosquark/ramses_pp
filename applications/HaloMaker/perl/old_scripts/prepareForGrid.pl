#!/usr/bin/perl
use strict;
use warnings;

my $fileRoot  = "/data/blaizot/TMP/Paulina/";
my $fileRoot2 = "/";

my $minFile   = 0;
my $maxFile   = 97;
my @snapList  = (0,1,2,9,11,14,20,25,33,34,37,39,45,48,55,57,58,65,69,73,76,79,97);
my $runDir    = "/data/blaizot/Horizon/1024-100h-1Mpc-W3/HaloMaker/FOF/";
my $paramFile = "/data/blaizot/Horizon/1024-100h-1Mpc-W3/HaloMaker/FOF/HaloMaker-svn/input_HaloMaker.dat";
my $execFile  = "/data/blaizot/Horizon/1024-100h-1Mpc-W3/HaloMaker/FOF/HaloMaker-svn/HaloMaker";
my $cmdFile   = "/data/blaizot/Horizon/1024-100h-1Mpc-W3/HaloMaker/FOF/cmdFile.sh";
open(CFILE,">$cmdFile");
#for (my $i = $minFile; $i <= $maxFile; $i++)
foreach (@snapList)
{
    my $i   = $_;
    my $num = "$i";
    if ($i < 100) {$num = "0"."$num";}
    if ($i < 10)  {$num = "0"."$num";}
    my $dir = "$runDir"."$i";
    my $cmd = "mkdir $dir";
    system($cmd);
    $cmd    = "cp $paramFile $dir";
    system($cmd);
    $cmd    = "cp $execFile $dir";
    system($cmd);
    my $file = "$dir"."/inputfiles_HaloMaker.dat";
    open(JFILE,">$file");
    print JFILE "\'$fileRoot"."$num"."$fileRoot2"."$num\'  Gd   1  $num \n";
    close(JFILE);
    print CFILE "cd $dir && ./HaloMaker > log.out \n";
}
close(CFILE);




