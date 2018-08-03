#!/usr/bin/perl

$LS    = 14;
$LT    = 14;
$H     = 10;
$m     = 0.00;
$beta  = 8.30;

$basedir = sprintf("/net/pool1/buividovich/toverlap/bin/b8.3000_s14_t14");
$datadir = sprintf("%s/H%i",$basedir, $H);

$wd      = sprintf("/net/pool1/buividovich/toverlap/bin/wd_su3_%i%id",   $LS, $LT);
$ovd     = sprintf("/net/pool1/buividovich/toverlap/bin/ovd_su3_%i%id", $LS, $LT);

print "\n\n\tDatadir: $datadir\n\n"; 

system("mkdir -p $datadir");

@configs = glob("$basedir/$LTx$LS*.dat");
$nconfigs = scalar(@configs);

print("$nconfigs configurations found in $basedir\n");

for($i=0; $i<$nconfigs; $i++)
{
 $config    = $configs[$i];
 $runfname  = sprintf("$basedir/run%03id_H%i.sh", $i, $H);
 $infofname = sprintf("$basedir/lwd_%03i.info", $i);
 $jobname   = sprintf("kalaydzhyan%03id", $i);
 $wdfname   = sprintf("$datadir/wd_%03id_s%i_t%i_H%i.dat", $i, $LS, $LT, $H);
 $ovdfname = sprintf("$datadir/ovd_%03id_%i.dat", $i, $H);
     
 
 open(INFOFILE,">$infofname");
 print INFOFILE "format    1\n";
 print INFOFILE "precision 4\n";
 print INFOFILE "lattice   $LS $LS $LS $LT\n";
 print INFOFILE "datafile  $config\n";
 print INFOFILE "sun       3\n";
 print INFOFILE "endian    l\n";
 close(INFOFILE);

 open(RUNFILE,">$runfname");
 print RUNFILE "\#PBS -N $jobname\n";
 print RUNFILE "\#PBS -q long\n";
 print RUNFILE "\#PBS -l cput=199:59:59,mem=900mb\n\n";
 print RUNFILE "cd $basedir\n";
# printf RUNFILE "$wd  -i $infofname -o $wdfname -E 0 -H $H -h -s -r 1.4 -a -S 30 -V\n";
 printf RUNFILE "$ovd -i $infofname -O $wdfname -o $ovdfname -E 0 -H $H -r 1.4 -m $m -t 0 -S 20 -V -a\n";
 close(RUNFILE);
 
 system("chmod a+x $runfname");
 
 system("qsub $runfname");
# system("renice 20 -u buividovich");
};

