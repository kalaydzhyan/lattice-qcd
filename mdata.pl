#!/usr/bin/perl

$beta          = 5.9;
$spacing       = 0.10;#
$LS            = 16;
$LT            = 4;
@Hs            = (   1,   2,   3,   5,  10,  15,   1,   2,   3,   4,   5,   6,   9,  10,  12,  15);
@Mus           = ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1);
$executable    = "/home/pools/1/buividovich/overlap/bin/mpl_su3_mu_164d";
$basedir       = "/home/pools/1/buividovich/overlap/bin/b5.9_s16_t4/";

$ns = scalar(@Hs);

for($i=0; $i<$ns; $i++)
{
 $H  = $Hs[$i];
 $Mu = $Mus[$i]; 
 $fsearch = sprintf("%s/H%i_mu%2.4f/ovd*.dat", $basedir, $H, $Mu); 
 $outfile = sprintf("mdata_b%2.4f_s%i_t%i_H%i_Mu%2.4f\n", $beta, $LS, $LT, $H, $Mu);
 @files = glob($fsearch);
 foreach $file (@files)
 {
  print "$file >> $outfile\n";
  system("$executable -i $file >> $outfile");
 };
};
