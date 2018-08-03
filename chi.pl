#!/usr/bin/perl

@sh     = ();
@s01    = ();
@ds01   = ();

open(S01, "<$ARGV[0]");
while(<S01>)
{
 @parts = split;
 push(@sh,   abs($parts[0]));
 push(@s01,  abs($parts[1]));
 push(@ds01, abs($parts[2]));
};
close(S01);

@sgmh   = ();
@sgm    = ();
@dsgm   = ();

open(SGM, "<$ARGV[1]");
while(<SGM>)
{
 @parts = split;
 push(@sgmh, abs($parts[0]));
 push(@sgm,  abs($parts[1]));
 push(@dsgm, abs($parts[2]));
};
close(SGM);

$nsh    = scalar(@sh);
$nsgmh  = scalar(@sgmh);

if($nsh!=$nsgmh)
{
 print "nsh ($nsh) != nsgmh ($nsgmh)!!!\n";
 exit;
};


open(CHI,">$ARGV[2]");
for($i=0; $i<$nsh; $i++)
{
 printf "At H=%2.4f GeV^2 Sigma = %2.4f MeV\n", $sgmh[$i], 1000*($sgm[$i]**0.333);
 $chi1 = ($s01[$i] + $ds01[$i])/($sgm[$i] - $dsgm[$i]);
 $chi2 = ($s01[$i] - $ds01[$i])/($sgm[$i] + $dsgm[$i]);
 $chi  = 0.5*($chi1 + $chi2);
 $dchi = 0.5*abs($chi1 - $chi2);
 printf "At H=%2.4f GeV^2 Chi = %2.4f +/- %2.4f\n", $sgmh[$i], $chi, $dchi;
 printf CHI "%2.4f %2.4f %2.4f\n", $sgmh[$i], $chi, $dchi;
};
close(CHI);