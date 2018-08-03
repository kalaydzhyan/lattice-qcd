#! /usr/bin/perl -w

$LS =  14; $LT = 14; $VOL = $LS*$LS*$LS*$LT;
$ND =   4; $NC = 2;
$BLOCK = $VOL*$ND*$NC + 1;
$NEV = 0;
$csp = 0.226044;

$rho  = 0;
$trho = 0;
$lambda = 0;

$evreader = "/home/pools/2/lena/overlap_14x14_ext/bin/read_ev_1414s";

$datadir = "/home/pools/2/lena/overlap_14x14_ext/bin/b3.2810_s14_t14/field_E=0_H=10/ovd_*.dat";
@files   = glob($datadir);
$nfiles  = scalar(@files);

@AX   = (0, 0, 0, 0);
@AX2  = (0, 0, 0, 0);
$maxlambda = 10;
$minlambda = 0;
$neiv = 0;

print "Total $nfiles files found.\n";

for($ifile = 0; $ifile<$nfiles; $ifile++)
{
 print "\n\n Processing the file $files[$ifile] ... \n\n";
 system("$evreader -i $files[$ifile] -o ./eigtemp.dat");
 open(TEMP,"<./eigtemp.dat");
 $counter = 0;
 while(<TEMP>)
 {
  @XY = split;
  if($counter%$BLOCK==0)
  {
   $dnorm = abs(sqrt(($XY[0] - 1)*($XY[0] - 1) + $XY[1]*$XY[1]) - 1);
   print "$lambda $trho\n";
   $trho  = 0;
   if($dnorm>0.00001)
   {
#    printf "Bad eigenvalue Nr. %i, dnorm = %2.4f: lies off the circle\n", $NEV, $dnorm;
#    printf "Re: %2.4f Im: %2.4f delta: %2.4f \n", $XY[0], $XY[1], $dnorm;
   };
   $lambda = abs(2*440*$XY[1]/($csp*(2.0 - $XY[0])) );
   if(($lambda >= $minlambda) && ($lambda < $maxlambda))
   {
    if($neiv>0)
    {
     print "\n\n ********************************************************* \n\n";
     print "With $neiv configurations: \n";
     for($i = 0; $i<4; $i++)
     {
      $DX = sqrt($AX2[$i]/$neiv - $AX[$i]*$AX[$i]/($neiv*$neiv));
      print "$i: $DX\n";
     };
     print "\n\n ********************************************************* \n\n";
    };
    $neiv++; 
   };
  }
  elsif(($lambda >= $minlambda) && ($lambda < $maxlambda))
  {
   $ixc = ($counter % $BLOCK) - 1;
   $idx = int($ixc/($ND*$NC)); 
   if($ixc % ($ND*$NC) == 0)
   {
    $trho += $rho;
    $tidx = $idx;
    @xc = (0, 0, 0, 0);
    for($i = 0; $i<4; $i++)
    {
     $xc[$i]   = $tidx % $LS;
     $tidx     = int($tidx/$LS);
     $AX[$i]  += $rho*$xc[$i];
     $AX2[$i] += $rho*$xc[$i]*$xc[$i];
    };
    $rho   = $XY[0]*$XY[0] + $XY[1]*$XY[1]; 
   }
   else
   {
    $rho  += $XY[0]*$XY[0] + $XY[1]*$XY[1];
   };
  };
  $counter++;
 };
 close(TEMP);
 system("rm -f ./bin/eigtemp.dat");
};

