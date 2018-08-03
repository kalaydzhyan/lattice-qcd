#!/usr/bin/perl -w

$LS   = 14;
$csp  = 0.208859;

$VOL3D = $LS*$LS*$LS;

$maxH = 10;

$mpl = "/home/pools/1/buividovich/overlap/bin/mpl_1414s";

@rho = ();

$maxrho = 0;

for($H=0; $H<=$maxH; $H++)
{
 $file = sprintf("/home/pools/1/buividovich/zmm/ovd_b3.2810_s14_t14_00001_sr_1.4_30_H$H.dat");
 print "Analyzing the file $file ... \n";
 
 open(READER, "$mpl -i $file -c $csp -n 6 -t 8|");
 while(<READER>)
 {
  push(@rho, $_);
  if($_ > $maxrho)
  {
   $maxrho = $_;
  };
 };
 close(READER);
};

$nrho = scalar(@rho);
$mrho = ($maxH + 1)*$VOL3D;
printf "$nrho $mrho %6.6f\n", $maxrho;

for($i=0; $i<$nrho; $i++)
{
 $rho[$i] = $rho[$i]/$maxrho;
};

for($H=0; $H<=$maxH; $H++)
{
 open(OUTFILE, ">sH$H.dat");
 for($i=0; $i<$VOL3D; $i++)
 {
  printf OUTFILE "%2.4f\n", $rho[$H*$VOL3D + $i];
 };
 close(OUTFILE);
 open(GENFNAME, ">sH$H.general");
 print GENFNAME "file = sH$H.dat\n";
 print GENFNAME "grid = $LS x $LS x $LS\n";
 print GENFNAME "format = ascii\n";
 print GENFNAME "interleaving = field\n";
 print GENFNAME "majority = row\n";
 print GENFNAME "field = field0\n";
 print GENFNAME "structure = scalar\n";
 print GENFNAME "type = float\n";
 print GENFNAME "dependency = positions\n";
 print GENFNAME "positions = regular, regular, regular, 0, 1, 0, 1, 0, 1\n";
 print GENFNAME "\n";
 print GENFNAME "end\n";
 close(GENFNAME);
};

system("zip -r movie.zip . -i sH*.dat -i sH*.general");
system("rm -f sH*.dat");
system("rm -f sH*.general");

