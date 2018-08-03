#!/usr/bin/perl

#for($hh=2; $hh<=7; $hh++)
#{

$LS    = 14;
$LT    = 14;
$H     = 1;
$E     = 0;
$m     = 0.00;
$beta  = 8.3;

$basedir  = sprintf("/afs/desy.de/user/t/tigrank/lattice/su3_over/bin/b%2.4f_s%i_t%i", $beta, $LS, $LT);
$datadir1 = sprintf("$basedir/..");
$scpdir   = sprintf("/net/pool1/buividovich/toverlap/bin/b%2.4f_s%i_t%i_single/H%i", $beta, $LS, $LT, $H);


$wd      = sprintf("./wd_su3_%i%is", $LS, $LT);
$ovd     = sprintf("./ovd_su3_%i%is", $LS, $LT);

print "\n\n\tDatadir: $datadir\n\n"; 

$basename = sprintf("%s/%ix%i*.dat",$basedir, $LT, $LS);

@configs = glob("$basename");
$nconfigs = scalar(@configs);

print("$nconfigs configurations found in $basedir\n");

for($i=0; $i<$nconfigs; $i++)
{
 
 $datadir2 = sprintf("%s/b%2.4f_s%i_t%i", $datadir1, $beta, $LS, $LT);
 $datadir  = sprintf("%s/H%i", $datadir2, $H);
   
 $config       = $configs[$i];
 $runfname     = sprintf("%s/run_ovd_%03is_H%i.sh", $basedir ,$i, $H);
 $infofname    = sprintf("%s/lwd_%03i.info", $basedir, $i);
 $jobname      = sprintf("kalaydzhyan%03is", $i);
 $wdfname      = sprintf("%s/wd_%03is_s%i_t%i_H%i.dat", $datadir1, $i, $LS, $LT, $H);
 $wdmarkname   = sprintf("wd_%03is_s%i_t%i_H%i", $i, $LS, $LT, $H);
 $ovdfname     = sprintf("%s/ovd_%03is.dat", $datadir, $i);
 
 open(INFOFILE,">$infofname");
 printf INFOFILE "format    1\n";
 printf INFOFILE "precision 4\n";
 printf INFOFILE "lattice   %i %i %i %i\n", $LS, $LS, $LS, $LT;
 printf INFOFILE "datafile  $config\n";
 printf INFOFILE "sun       3\n";
 printf INFOFILE "endian    l\n";
 close(INFOFILE);

 open(RUNFILE,">$runfname");
 printf RUNFILE "\#\$ -l h_rt=125:00:00\n";
 printf RUNFILE "\#\$ -l h_vmem=1G\n";
 printf RUNFILE "\#\$ -l arch=x86\n\n";
 
 
 printf RUNFILE "mkdir -p $datadir1\n";
 printf RUNFILE "mkdir -p $datadir2\n";
 printf RUNFILE "mkdir -p $datadir\n";
 
 printf RUNFILE "cd $basedir\n";
 printf RUNFILE "cd ..\n";

# printf RUNFILE "while [ -f $basedir/../wait ]\n";
# printf RUNFILE "do\n";
# printf RUNFILE "sleep 30\n";
# printf RUNFILE "done\n";
# printf RUNFILE "cp $basedir/../mark $basedir/../wait\n";
# printf RUNFILE "cd ..\n";
# printf RUNFILE "./compile.sh %i %i\n", $LS, $LT;
# printf RUNFILE "cd bin\n";
# printf RUNFILE "rm -f $basedir/../wait\n";

# printf RUNFILE "$wd  -i $infofname -o $wdfname -E $E -H $H -h -s -r 1.4 -a -S 30 -V\n";

 printf RUNFILE "while [ -f $basedir/../wait ]\n";
 printf RUNFILE "do\n";
 printf RUNFILE "sleep 30\n";
 printf RUNFILE "done\n";
 printf RUNFILE "cp $basedir/../mark $basedir/../wait\n";
 printf RUNFILE "cd ..\n";
 printf RUNFILE "./compile.sh %i %i 1\n", $LS, $LT;
 printf RUNFILE "scp buividovich\@rrcsrv.itep.ru:$scpdir/wd_%03is* $basedir/../\n", $i;
 printf RUNFILE "cd bin\n";
 printf RUNFILE "cp $basedir/../mark $basedir/../$wdmarkname\n";
 printf RUNFILE "rm -f $basedir/../wait\n";

 printf RUNFILE "$ovd -i $infofname -O $wdfname -o $ovdfname -E $E -H $H -r 1.4 -m $m -t 0 -S 20 -V -a\n";
 printf RUNFILE "scp $ovdfname buividovich\@rrcsrv.itep.ru:$scpdir/\n";
 printf RUNFILE "rm -f $wdfname\n";
 printf RUNFILE "rm -f $ovdfname\n";
 close(RUNFILE);
 
 system("chmod a+x $runfname");

 $markname = sprintf("%s/../wait*",$basedir);
 @marks = glob("$markname");
 $nmarks = scalar(@marks);
# system("echo $nmarks");
 
 while($nmarks != 0)
 {
  $markname = sprintf("%s/../wait*",$basedir);
  @marks = glob("$markname");
  $nmarks = scalar(@marks);
 sleep(5);
 }
 
system("qsub -cwd $runfname");
sleep(15);
};

#}