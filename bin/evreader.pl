#! /usr/bin/perl -w

$LS =  4; $LT = 4;
$ND = 4;  $NC = 2;

#$CSP = 0.226044; # - for 14x14, beta = 3.2810
$CSP = 0.397594; # - for 4x4, beta = 2.95
#$CSP  =1; # for testing purposes 

#$evreader = "/home/pools/2/lena/overlap_14x14_ext/bin/read_ev_1414s";
$evreader = "/home/pools/1/buividovich/overlap/bin/read_ev_44s";

$ARGC = scalar(@ARGV);
die "correct usage: evreader.pl file [outfile] \n extracts eigenvalues from data files produced by overlap program and appends them to the outfile, if the latter is specified \n" unless ($ARGC==1 || $ARGC==2);

$VOL = $LS*$LS*$LS*$LT;
$BLOCK = $VOL*$ND*$NC + 1;

system("$evreader -i $ARGV[0] -o ./eigtemp.dat");
open(TEMP,"<./eigtemp.dat");

$counter   = 0;
$evcounter = 0;

if($ARGC==2)
{
 open(OUTFILE, ">>$ARGV[1]");
};

while(<TEMP>)
{
 if($counter%$BLOCK==0)
 {
  @XY = split;
# $dnorm = abs(sqrt(($XY[0] - 1)*($XY[0] - 1) + $XY[1]*$XY[1]) - 1);
# if($dnorm>0.00001)
# {
#  printf "Bad eigenvalue Nr. %i, dnorm = %2.4f: lies off the circle\n", $evcounter+1, $dnorm;
#  printf "Re: %2.4f Im: %2.4f delta: %2.4f \n", $XY[0], $XY[1], $dnorm;
# };
# $lambda = abs( 2*440.0*$XY[1]/($CSP*(2.0 - $XY[0])) );
# printf "Eigenvalue No %i: %4.4f\n", $evcounter, $lambda;
# push(@EVS2,$lambda);
  $evcounter++;
  if($ARGC==2)
  {
   printf OUTFILE "%8.8f %8.8f\n", $XY[0], $XY[1];
  }
  else
  {
   printf "%4.6f %4.6f\n", $XY[0], $XY[1];
  };
# $NEV++;
 };
 $counter++; 
};
close(TEMP);
system("rm -f ./eigtemp.dat");
print "Total $evcounter eigenvalues read\n\n";

if($ARGC==2)
{
 close(OUTFILE);
};

#@EVS2 = sort {$a <=> $b} @EVS2;

#@EVS = ();
#for($i = 0; $i<($NEV/2); $i++)
#{
# $diff = abs($EVS2[2*$i] - $EVS2[2*$i + 1]);
# if($diff>0.0001)
# {
#  printf "Bad eigenvalues Nr. %i, diff = %2.4f: unpaired\n", $NEV, $dnorm;
# };
# push(@EVS,$EVS2[2*$i]);
#};

#print "Sorted eigenvalues:\n\n";
#for($i = 0; $i<($NEV/2); $i++)
#{
# printf "%02i: %4.4f MeV\n", ($i+1), $EVS[$i]; 
#};
