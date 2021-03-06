#!/usr/bin/perl

$maxnfiles     = 16;

$beta          = 3.3555;
$spacing       = 0.089;
$mass_phys     = 50.0; # IR Cutoff - in MeV
$max_nz_evals  = 10;    # UV cutoff
$LS            = 16;
$LT            = 16;
#@Hs	       = (0);
@Hs	       = (0, 5, 10);
$executable    = "./mcurr_1616s"; 
$basedir       = "/net/pool1/buividovich/overlap/bin/b3.3555_s16_t16";

system("renice 20 -u buividovich");

$NH = scalar(@Hs);
for($i=0; $i<$NH; $i++)
{
 $H = $Hs[$i];

 $paramname = sprintf("params_%1.4f_s%i_t%i_H%i.in", $beta, $LS, $LT, $H);

 open(PARAMS, ">$paramname");
 printf PARAMS "%2.4f\n",  $beta;
 printf PARAMS "%2.4f\n",  $spacing;
 printf PARAMS "%i\n",     $H;
 printf PARAMS "%4.4lf\n", $mass_phys;
 printf PARAMS "%i\n",     $max_nz_evals;
 printf PARAMS "b%2.4f_H%i_s%i_t%i_m%i_ev%i\n", $beta, $H, $LS, $LT, int($mass_phys), $max_nz_evals;
 close(PARAMS);
 
 @files = glob("$basedir/H$H/ovd*.dat");
 
 $nfiles = scalar(@files);
 $nfiles = ($nfiles > $maxnfiles)? $maxnfiles : $nfiles;
 
 @files = @files[0..$nfiles];
 $arg = join(" ", @files);

 $scriptname = sprintf("shedule_%1.4f_s%i_t%i_H%i.sh", $beta, $LS, $LT, $H);

 system("rm $scriptname");
 open(TASK, ">>$scriptname");
 printf TASK "cd /home/pools/1/buividovich/toverlap/bin/\n";
 printf TASK "cp $paramname params.in\n";
 printf TASK "$executable $arg\n";
 printf TASK "rm $paramname";
 close(TASK);
 system("qsub $scriptname");
 print "batch with H=$H has been added\n";

};

#system("rm shedule");