#! /usr/bin/perl -w

$ARGC = scalar(@ARGV);
die "correct usage: current_plotter.pl evn infile outfile\n extracts configuration of currents (color indices contracted) of eigenmode number evn from data files produced by overlap program \n" unless $ARGC==3;

$evreader = "/home/pools/2/lena/overlap_14x14_ext/bin/read_ev_1414s";

$LS =  14; $LT = 14; $VOL = $LS*$LS*$LS*$LT;
$ND =   4; $NC = 2;
$BLOCK = $VOL*$ND*$NC + 1;
$NEV = 0;
$csp = 0.226044;

$sps = $LS*$LS*$LS;     # Sites per slice
$jsc = 1000000;           # A coefficient to rescale the currents

$j1 = 0; $j2 = 0; $j3 = 0; # Currents. J0 is not considered, since it is equal to mode density

@psi_r = (0, 0, 0, 0); # real and imaginary components of the wave function 
@psi_i = (0, 0, 0, 0);

print "\n\n Processing the file $ARGV[1] ... \n\n";
system("$evreader -i $ARGV[1] -o ./eigtemp.dat");
system("mkdir -p ./jslices");

open(TEMP,"<./eigtemp.dat");

$counter   = 0;
$evcounter = 0;

while(<TEMP>)
{
 @XY = split;
 if($counter%$BLOCK==0)
 {
  $dnorm = abs(sqrt(($XY[0] - 1)*($XY[0] - 1) + $XY[1]*$XY[1]) - 1);
  if($evcounter == $ARGV[0])
  {
   $j1 += 2*($psi_i[0]*$psi_i[1] - $psi_i[2]*$psi_i[3] + $psi_r[0]*$psi_r[1] - $psi_r[2]*$psi_r[3]);
   $j2 += 2*($psi_i[1]*$psi_r[0] - $psi_i[0]*$psi_r[1] - $psi_i[3]*$psi_r[2] + $psi_i[2]*$psi_r[3]);
   $j3 += $psi_i[0]*$psi_i[0] - $psi_i[1]*$psi_i[1] - $psi_i[2]*$psi_i[2] + $psi_i[3]*$psi_i[3] + $psi_r[0]*$psi_r[0] - $psi_r[1]*$psi_r[1] - $psi_r[2]*$psi_r[2] + $psi_r[3]*$psi_r[3];  
   $j1 *= $jsc;
   $j2 *= $jsc;   
   $j3 *= $jsc;
   printf OUTDATA "%4.2f %4.2f %4.2f\n", $j1, $j2, $j3;
  };
  print "$evcounter: $lambda\n";
  if($dnorm>0.00001)
  {
    printf "Bad eigenvalue Nr. %i, dnorm = %2.4f: lies off the circle\n", $NEV, $dnorm;
    printf "Re: %2.4f Im: %2.4f delta: %2.4f \n", $XY[0], $XY[1], $dnorm;
  };
  $lambda = abs(2*440*$XY[1]/($csp*(2.0 - $XY[0])) );
  $evcounter++;
  if($evcounter == $ARGV[0])
  {
   printf "Analyzing eigenmode No %i, eigenvalue %4.6f ... \n", $evcounter, $lambda;
  };
 } 
 elsif($evcounter == $ARGV[0])
 {
  $ixc = ($counter % $BLOCK) - 1;
  $idx = int($ixc/($ND*$NC)); 
  $dix = $ixc % $ND; # Dirac spinor index
  if($dix == 0)
  {
   $j1 += 2*($psi_i[0]*$psi_i[1] - $psi_i[2]*$psi_i[3] + $psi_r[0]*$psi_r[1] - $psi_r[2]*$psi_r[3]);
   $j2 += 2*($psi_i[1]*$psi_r[0] - $psi_i[0]*$psi_r[1] - $psi_i[3]*$psi_r[2] + $psi_i[2]*$psi_r[3]);
   $j3 += $psi_i[0]*$psi_i[0] - $psi_i[1]*$psi_i[1] - $psi_i[2]*$psi_i[2] + $psi_i[3]*$psi_i[3] + $psi_r[0]*$psi_r[0] - $psi_r[1]*$psi_r[1] - $psi_r[2]*$psi_r[2] + $psi_r[3]*$psi_r[3];  
  };
  if($ixc % ($ND*$NC) == 0)
  {
#********************************************************#
# Here we have the value of mode density at the site idx #
#********************************************************#
   if($ixc!=0)
   {
    $j1 *= $jsc;
    $j2 *= $jsc;
    $j3 *= $jsc;
    printf OUTDATA "%4.2f %4.2f %4.2f\n", $j1, $j2, $j3;
    $j1 = 0; $j2 = 0; $j3 = 0;
   };
   if($idx%$sps == 0)
   {
    $slice = int($idx/$sps);
    $outfname = "./jslices/js$slice.dat";
    $genfname  ="./jslices/js$slice.general";
    print "outfname: $outfname, slice: $slice, idx: $idx, LS: $LS\n";
    if($ixc!=0)
    {
     close(OUTDATA);
    };
    open(GENFNAME, ">$genfname");
    print GENFNAME "file = js$slice.dat\n";
    print GENFNAME "grid = $LS x $LS x $LS\n";
    print GENFNAME "format = ascii\n";
    print GENFNAME "interleaving = field\n";
    print GENFNAME "majority = row\n";
    print GENFNAME "field = field0\n";
    print GENFNAME "structure = 3-vector\n";
    print GENFNAME "type = float\n";
    print GENFNAME "dependency = positions\n";
    print GENFNAME "positions = regular, regular, regular, 0, 1, 0, 1, 0, 1\n";
    print GENFNAME "\n";
    print GENFNAME "end\n";
    close(GENFNAME);
    open(OUTDATA,">$outfname");
   };
  };
  $psi_r[$dix] = $XY[0];
  $psi_i[$dix] = $XY[1];
 };
 $counter++;
};
close(OUTDATA);

system("rm -f ./eigtemp.dat");
system("zip -r -j $ARGV[2] ./jslices");
system("rm -r -f ./jslices");
