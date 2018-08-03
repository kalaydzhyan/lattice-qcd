#!/usr/bin/perl

$maxnfiles     = 30;

$beta          = 3.2810;
$spacing       = 0.102747;
$mass_phys     = 50.0; # IR Cutoff - in MeV
$max_nz_evals  = 12;    # UV cutoff
$LS            = 14;
$LT            = 14;
@Hs	       = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 20);
$executable    = "/home/pools/1/buividovich/overlap/bin/jxjy_1414s";
$basedir       = "/home/pools/1/buividovich/overlap/bin/b3.2810_s14_t14/";

$NH = scalar(@Hs);
for($i=0; $i<$NH; $i++)
{
 $H = $Hs[$i];
 open(PARAMS, ">params.in");
 printf PARAMS "%2.4f\n",  $beta;
 printf PARAMS "%2.4f\n",  $spacing;
 printf PARAMS "%i\n",     $H;
 printf PARAMS "%4.4lf\n", $mass_phys;
 printf PARAMS "%i\n",     $max_nz_evals;
 printf PARAMS "b%2.4f_mu_s%i_t%i_m%i_ev%i\n", $beta, $LS, $LT, int($mass_phys), $max_nz_evals;
 close(PARAMS);
 
 @files = glob("$basedir/H$H/ovd*.dat");
 
 $nfiles = scalar(@files);
 print "$nfiles files found, ";
 $nfiles = ($nfiles > $maxnfiles)? $maxnfiles : $nfiles;
 print "processing only the first $nfiles.\n";

 @files = @files[0..$nfiles];
 $arg = join(" ", @files);

 print "$arg\n";

 system("$executable $arg");
};

