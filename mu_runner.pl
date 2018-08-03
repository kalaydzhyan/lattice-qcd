#!/usr/bin/perl

$mu = 0.0;

$execpath = "/home/pools/1/buividovich/overlap/bin/";
$datapath = "/home/pools/1/buividovich/overlap/bin/b2.95_s4_t4/";
$outpath  = "/home/pools/1/buividovich/overlap/";

$wdname   = $execpath."wd_44s";
$ovdname  = $execpath."ovd_44s";
$evreader = $execpath."evreader.pl";

$outfile  = $outpath."eigvals_mu00.dat";

for($ifile = 1; $ifile<=9; $ifile++)
{
 open(RUNFILE, ">run.sh");
 
 $infile     = sprintf("%slwd_%05i.info", $datapath, $ifile);
 $wdoutfile  = sprintf("%swd_b2.9500_s4_t4_%05i_sr_1.4_30.dat", $datapath, $ifile);
 $ovdoutfile = sprintf("%sovd_b2.9500_s4_t4_%05i_sr_1.4_30.dat", $datapath, $ifile);
 
 print "infile: $infile\n";
 print "wdoutfile: $wdoutfile\n";
 print "ovdoutfile: $ovdoutfile\n"; 
 
 print "\n\n";
 
 printf RUNFILE "cd %s\n", $datapath;
 printf RUNFILE "%s -i %s -o %s -M %2.4f -h -s -r 1.4 -S 30 -L 15 -a -F 5 -V -Q -R\n", $wdname, $infile, $wdoutfile, $mu;
 printf RUNFILE "%s -i %s -O %s -o %s -M %2.4f -r 1.4 -S 30 -t 0 -a -P 50 -L 0.0000001 -T 0.00001\n", $ovdname, $infile, $wdoutfile, $ovdoutfile, $mu;

# the commands below are for the Morozov's version of the code 
# printf RUNFILE "%s -i %s -o %s -h -s -r 1.4 -S 30 -L 2 -a -F 3\n", $wdname, $infile, $wdoutfile;
# printf RUNFILE "%s -i %s -O %s -o %s -r 1.4 -S 30 -t 0 -a\n", $ovdname, $infile, $wdoutfile, $ovdoutfile;

 printf RUNFILE "%s %s %s\n", $evreader, $ovdoutfile, $outfile;
 
 printf RUNFILE "rm -f $wdoutfile\n";
 printf RUNFILE "rm -f $ovdoutfile\n";
 printf RUNFILE "rm -f %sarnoldi.log\n", $execpath;
 printf RUNFILE "rm -f %sarnoldi.log\n", $datapath;
 
 close(RUNFILE);
 
 system("chmod a+x ./run.sh");
 system("./run.sh");
 sleep(60);
};
