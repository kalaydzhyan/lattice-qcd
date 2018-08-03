#! /usr/bin/perl -w

$ARGC = scalar(@ARGV);
die "correct usage: mode_plotter.pl evn infile outfile\n extracts eigenmode number evn from data files produced by overlap program \n" unless $ARGC==3;

#$evreader = "/home/pools/2/lena/overlap_14x14_ext/bin/read_ev_1414s";
$evreader = "/home/pools/1/buividovich/overlap/bin/read_ev_66s";

$LS =  6; $LT = 6; $VOL = $LS*$LS*$LS*$LT;
$ND =   4; $NC = 2;
$BLOCK = $VOL*$ND*$NC + 1;
$NEV = 0;
#$csp = 0.226044;
$csp = 1;

$sps = $LS*$LS*$LS;     # Sites per slice
$rho  = 0;
$ipr   = 0;

print "\n\n Processing the file $ARGV[1] ... \n\n";
system("$evreader -i $ARGV[1] -o ./eigtemp.dat");
system("mkdir -p ./slices");

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
   $trho += $rho;
   $ipr  += $rho*$rho;
   $grho = 1000000*$rho;
   printf OUTDATA "%2.4f\n", $grho;
  };
  print "$lambda $trho\n";
  $trho  = 0;
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
  if($ixc % ($ND*$NC) == 0)
  {
#********************************************************#
# Here we have the value of mode density at the site idx #
#********************************************************#
   $trho += $rho;
   $ipr  += $rho*$rho;
   $grho = 1000000*$rho;
   if($ixc!=0)
   {
    printf OUTDATA "%4.2f\n", $grho;
   };
   if($idx%$sps == 0)
   {
    $slice = int($idx/$sps);
    $outfname = "./slices/s$slice.dat";
    $genfname  ="./slices/s$slice.general";
    print "outfname: $outfname, slice: $slice, idx: $idx, LS: $LS\n";
    if($ixc!=0)
    {
     close(OUTDATA);
    };
    open(GENFNAME, ">$genfname");
    print GENFNAME "file = s$slice.dat\n";
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
    open(OUTDATA,">$outfname");
   };
#********************************************************#   
# Here we start calculating rho again                    #
#********************************************************#
   $rho   = $XY[0]*$XY[0] + $XY[1]*$XY[1]; 
  }
  else
  {
   $rho  += $XY[0]*$XY[0] + $XY[1]*$XY[1];
  };
 };
 $counter++;
};
close(OUTDATA);

$ipr = $VOL*$ipr;
printf "For the eigenmode under analysis, the IPR is: %6.6f\n", $ipr;

system("rm -f ./eigtemp.dat");
system("zip -r -j $ARGV[2] ./slices");
system("rm -r -f ./slices");
