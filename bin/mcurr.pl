#!/usr/bin/perl

$maxnfiles     = 30;
$usemu	       = 0; #1 if chemical potential included, zero otherwise

$beta          = 3.2810;
$spacing       = 0.102747;
$mass_phys     = 50.0; # IR Cutoff - in MeV
$max_nz_evals  = 12;    # UV cutoff
$LS            = 14;
$LT            = 14;
@Hs	       = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 20);
$executable    = "/home/pools/1/buividovich/toverlap/bin/check_1414s";
$basedir       = "/home/pools/1/buividovich/overlap/bin/b3.2810_s14_t14/";

#$beta          = 5.9;
#$spacing       = 0.10;#
#$mass_phys     = 50.0; # IR Cutoff - in MeV
#$max_nz_evals  = 7;    # UV cutoff
#$LS            = 16;
#$LT            = 4;
#@Hs	       = (   1,   2,   3,   5,  10,  15,   1,   2,   3,   4,   5,   6,   9,  10,  12,  15);
#@Mus	       = ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1);
#$executable    = "/home/pools/1/buividovich/overlap/bin/mcurr_su3_mu_164d";
#$basedir       = "/home/pools/1/buividovich/overlap/bin/b5.9_s16_t4/";

#$beta          = 3.3250;
#$spacing       = 0.0949361;
#$mass_phys     = 50.0; # IR Cutoff - in MeV
#$max_nz_evals  = 10;    # UV cutoff
#$LS            = 16;
#$LT            = 6;
#@Hs	       = (0, 1, 2, 3, 5, 10);
#$executable    = "/home/pools/1/buividovich/overlap/bin/mcurr_166s";
#$basedir       = "/home/pools/1/buividovich/overlap/bin/b3.3250_s16_t6/";

#$beta          = 3.1600;
#$spacing       = 0.126044;
#$mass_phys     = 50.0; # IR Cutoff - in MeV
#$max_nz_evals  = 10;    # UV cutoff
#$LS            = 16;
#$LT            = 6;
#@Hs	       = (0, 1, 2, 3, 5, 8, 10, 12);
#@Hs	       = (5);
#$executable    = "/home/pools/1/buividovich/overlap/bin/mcurr_166s";
#$executable    = "/home/pools/1/buividovich/overlap/bin/cr52j2_166s";
#$basedir       = "/home/pools/1/buividovich/overlap/bin/b3.1600_s16_t6/";

#$beta          = 3.3555;
#$spacing       = 0.0894159;
#$mass_phys     = 50.0; # IR Cutoff - in MeV
#$max_nz_evals  = 12;    # UV cutoff
#$LS            = 16;
#$LT            = 16;
#@Hs	       = (0, 5, 10);
#$executable    = "/home/pools/1/buividovich/overlap/bin/mcurr_1616s";
#$basedir       = "/home/pools/1/buividovich/overlap/bin/b3.3555_s16_t16/";

#$beta          = 3.2810;
#$spacing       = 0.102747;
#$mass_phys     = 50.0; # IR Cutoff - in MeV
#$max_nz_evals  = 12;    # UV cutoff
#$LS            = 16;
#$LT            = 16;
#@Hs	        = (0, 3, 5);
#$executable    = "/home/pools/1/buividovich/overlap/bin/mcurr_1616s";
#$basedir       = "/home/pools/1/buividovich/overlap/bin/b3.2810_s16_t16/";

#$beta          = 3.2810;
#$spacing       = 0.102747;
#$mass_phys     = 50.0; # IR Cutoff - in MeV
#$max_nz_evals  = 10;    # UV cutoff
#$LS            = 8;
#$LT            = 8;
#@Hs	       = (0, 1, 2, 3, 4, 5);
#$executable    = "/home/pools/1/buividovich/overlap/bin/mcurr_88s";
#$executable    = "/home/pools/1/buividovich/overlap/bin/cr52j2_88s";
#$basedir       = "/home/pools/1/buividovich/overlap/bin/instanton/";

#$beta          = 5.9;
#$spacing       = 0.102747;
#$mass_phys     = 50.0; # IR Cutoff - in MeV
#$max_nz_evals  = 7;    # UV cutoff
#$LS            = 8;
#$LT            = 8;
#@Hs	       = (0, 5);
#@Mus	       = (0.1, 0.1);	
#$executable    = "/home/pools/1/buividovich/overlap/bin/mcurr_su3_mu_44d";
#executable    = "/home/pools/1/buividovich/overlap/bin/cr52j2_88s";
#$basedir       = "/home/pools/1/buividovich/overlap/bin/b5.9_s4_t4/";

$NH = scalar(@Hs);
for($i=0; $i<$NH; $i++)
{
 $H = $Hs[$i];
 if($usemu==1)
 {
  $Mu = $Mus[$i];
 };
 open(PARAMS, ">params.in");
 printf PARAMS "%2.4f\n",  $beta;
 printf PARAMS "%2.4f\n",  $spacing;
 printf PARAMS "%i\n",     $H;
 if($usemu==1)
 {
  printf PARAMS "%2.4f\n", $Mu;
 };
 printf PARAMS "%4.4lf\n", $mass_phys;
 printf PARAMS "%i\n",     $max_nz_evals;
 if($usemu==0)
 {
  printf PARAMS "b%2.4f_s%i_t%i_m%i_ev%i\n", $beta, $LS, $LT, int($mass_phys), $max_nz_evals;
 }
 else
 {
  printf PARAMS "b%2.4f_mu_s%i_t%i_m%i_ev%i\n", $beta, $LS, $LT, int($mass_phys), $max_nz_evals;
 }; 
 close(PARAMS);
 
 if($usemu==1)
 {
  $fsearch = sprintf("%s/H%i_mu%2.4f/ovd*.dat", $basedir, $H, $Mu);
  @files = glob($fsearch);
 }
 else
 {
  @files = glob("$basedir/H$H/ovd*.dat");
 }; 

 $nfiles = scalar(@files);
 print "$nfiles files found, ";
 $nfiles = ($nfiles > $maxnfiles)? $maxnfiles : $nfiles;
 print "processing only the first $nfiles.\n";

 @files = @files[0..$nfiles];
 $arg = join(" ", @files);

 print "$arg\n";

 system("$executable $arg");
};

