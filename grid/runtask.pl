#!/usr/bin/perl

$LS    = 14;
$LT    = 14;
$H     = 1;
$E     = 0;
$m     = 0.00;
$beta  = 8.30;
$LOGNAME = sprintf("kalaydzhyan");

$basedir = sprintf("b%2.4f_s%i_t%i",$beta, $LS, $LT);
$datadir = sprintf("%s/H%i",$basedir, $H);

$wd      = sprintf("wd_su3_%i%id",  $LS, $LT);
$ovd     = sprintf("ovd_su3_%i%id", $LS, $LT);

print "\n\n\tDatadir: $datadir\n\n"; 

system("mkdir -p $datadir");
$configname = sprintf("%s/%ix%i*.dat", $basedir, $LT, $LS);

@configs = glob("$configname");
$nconfigs = scalar(@configs);

print("$nconfigs configurations found in $basedir\n");

for($i=0; $i<$nconfigs; $i++)
{
 $config      = $configs[$i];
 @configarray = split('/', "$config");
 $shortconfig = $configarray[scalar(@configarray)-1];
 $runfname    = sprintf("run%03id_H%i.sh", $i, $H);
 $jdlfname    = sprintf("task%03id_H%i.jdl", $i, $H);
 $infofname   = sprintf("lwd_%03i.info", $i);
 $wdfname     = sprintf("$datadir/wd_%03id_s%i_t%i_H%i.dat", $i, $LS, $LT, $H);
 $ovdfname    = sprintf("$datadir/ovd_%03id_s%i_t%i_H%i.dat", $i, $LS, $LT, $H);
 $ovdshortfname = sprintf("ovd_%03id_s%i_t%i_H%i.dat", $i, $LS, $LT, $H);
 
# system("lcg-cp -D srmv2 -b -v --vo lattice.itep.ru file:$config srm://selattice.itep.ru:8446/srm/managerv2?SFN=/dpm/itep.ru/home/lattice.itep.ru/$LOGNAME/$shortconfig");
 
 open(INFOFILE,">$infofname");
 print INFOFILE "format    1\n";
 print INFOFILE "precision 4\n";
 print INFOFILE "lattice   $LS $LS $LS $LT\n";
 print INFOFILE "datafile  $basedir/$shortconfig\n";
 print INFOFILE "sun       3\n";
 print INFOFILE "endian    l\n";
 close(INFOFILE);

 open(RUNFILE,">$runfname");
 print RUNFILE "\#PBS -q long\n";
 print RUNFILE "\#PBS -l cput=199:59:59,mem=900mb\n\n";
 print RUNFILE "tar -xvf toverlap.tar\n";
 print RUNFILE "cd toverlap\n";
 print RUNFILE "./compile.sh $LS $LT 0\n";
 print RUNFILE "cd bin\n";
 print RUNFILE "mkdir $basedir\n";
 print RUNFILE "mkdir $datadir\n";
 print RUNFILE "mv ../../$infofname ./$basedir\n";
 print RUNFILE "lcg-cp -D srmv2 -b -v --vo lattice.itep.ru srm://selattice.itep.ru:8446/srm/managerv2?SFN=/dpm/itep.ru/home/lattice.itep.ru/$LOGNAME/$shortconfig file:\$PWD/$basedir/$shortconfig\n";
 print RUNFILE "./$wd  -i $basedir/$infofname -o $wdfname -E $E -H $H -h -s -r 1.4 -a -S 30 -V\n";
 print RUNFILE "./$ovd -i $basedir/$infofname -O $wdfname -o $ovdfname -E $E -H $H -r 1.4 -m $m -t 0 -S 20 -V -a\n";
 print RUNFILE "lcg-cp -D srmv2 -b -v --vo lattice.itep.ru file:\$PWD/$ovdfname srm://selattice.itep.ru:8446/srm/managerv2?SFN=/dpm/itep.ru/home/lattice.itep.ru/$LOGNAME/$ovdshortfname\n";
 close(RUNFILE);
 
 open(JDLFILE, ">$jdlfname");
 print JDLFILE "VirtualOrganisation = \"lattice.itep.ru\";\n";
 print JDLFILE "Executable = \"$runfname\";\n";
 print JDLFILE "Arguments = \"\";\n";
 print JDLFILE "StdOutput = \"std.out\";\n";
 print JDLFILE "StdError = \"std.err\";\n";
 print JDLFILE "OutputSandbox = \{\"std.out\",\"std.err\"\};\n";
 print JDLFILE "InputSandbox = \{\"$runfname\",\"$infofname\", \"toverlap.tar\"\};\n";
 print JDLFILE "Requirements=RegExp\(\"ce3.itep.ru\",other.GlueCEUniqueID);\n";
 print JDLFILE "JobType = \"Normal\";\n";
 print JDLFILE "RetryCount = 7;\n";
 close(JDLFILE);
 
 system("chmod 777 *.sh");
 system("glite-wms-job-submit -o job_id -a $jdlfname");
 system("rm -f $jdlfname");
 system("rm -f $runfname");
 system("rm -f $infofname"); 
};

