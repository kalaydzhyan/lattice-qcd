#!/usr/bin/perl

do "/home/pools/1/buividovich/overlap/bin/mesons/runparams.pl";

$wd      = sprintf("/home/pools/1/buividovich/overlap/bin/wd_%i%id",   $LS, $LT);
$prop    = sprintf("/home/pools/1/buividovich/overlap/bin/prop_%i%id", $LS, $LT);

print "\n\n\tDatadir: $datadir\n\n"; 

system("mkdir -p $datadir");

@configs = glob("$basedir/CON*");
$nconfigs = scalar(@configs);

for($i=0; $i<$nconfigs; $i++)
{
 $config    = $configs[$i];
 $runfname  = sprintf("$basedir/run%03id.sh", $i);
 $infofname = sprintf("$basedir/lwd_%03i.info", $i);
 $jobname   = sprintf("cme%03id", $i);
 $wdfname   = sprintf("$datadir/wd_%03id.dat", $i);
 $propfname = sprintf("$datadir/prop_%03id.dat", $i);
 
 open(INFOFILE,">$infofname");
 print INFOFILE "format    2\n";
 print INFOFILE "precision 4\n";
 print INFOFILE "lattice   $LS $LS $LS $LT\n";
 print INFOFILE "datafile  $config\n";
 print INFOFILE "sun       2\n";
 print INFOFILE "endian    l\n";
 close(INFOFILE);

 open(RUNFILE,">$runfname");
 print RUNFILE "\#PBS -N $jobname\n";
 print RUNFILE "\#PBS -q long\n";
 print RUNFILE "\#PBS -l cput=199:59:59,mem=900mb\n\n";
 print RUNFILE "cd $basedir\n";
#printf RUNFILE "$wd -i $infofname -o $wdfname -E 0 -H $H -h -s -r 1.4 -a -S 30 -V\n";
 printf RUNFILE "$prop -i $infofname -O $wdfname -o $propfname -E 0 -H $H -r 1.4 -m $m -t 0 -a\n";
 close(RUNFILE);
 
 system("chmod a+x $runfname");
 
 system("qsub $runfname");
};

