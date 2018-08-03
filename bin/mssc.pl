#!/usr/bin/perl -w

$mssc = "/home/pools/1/buividovich/overlap/bin/mssc_166s";
$csp  = 0.208859;

die "Correct usage: mssc.pl indir\n" unless (scalar(@ARGV)==1); 

$indir = $ARGV[0];

@flist = glob "$indir/ovd*.dat";
$nfiles = scalar(@flist);

$m01  = 0.0;
$m01d = 0.0;
$m23  = 0.0;
$m23d = 0.0;
$j    = 0.0;
$jd   = 0.0;
$j5   = 0.0;
$j5d  = 0.0;
$nm   =   0;

$maxnfiles = ($nfiles > 20)? 20 : $nfiles;

print "\n\n\t $nfiles files found, processing the first $maxnfiles only.\n\n";

for($ifile=0; $ifile<$maxnfiles; $ifile++)
{
 $file = $flist[$ifile];
 print "Analyzing the file $file ... \n";
 
 $cm01  = 0.0;
 $cm01d = 0.0;
 $cm23  = 0.0;
 $cm23d = 0.0;
 $cj    = 0.0;
 $cjd   = 0.0;
 $cj5   = 0.0;
 $cj5d  = 0.0;
 $cnm   = 0;

 open(READER, "$mssc -i $file -c $csp|");
 
 $count  = 0;
 $zmodes = 0;
 while(<READER>)
 {
  @mps = split;
  die "Wrong output of mssc!!!\n" unless (scalar(@mps)==6);
  
  if($mps[0]<1 && $mps[5]>0.8)
  {
   printf "Zero mode! Lambda = %4.6f, chirality = %2.2lf, M01 = %4.6lf, M23 = %4.6lf, J = %4.6lf, J5 = %4.6lf\n", $mps[0], $mps[5], $mps[1], $mps[2], $mps[3], $mps[4];
   $cm01  += $mps[1];
   $cm23  += $mps[2];
   $cj    += $mps[3];
   $cj5   += $mps[4];
   $cm01d += $mps[1]*$mps[1];
   $cm23d += $mps[2]*$mps[2];
   $cjd   += $mps[3]*$mps[3];
   $cj5d  += $mps[4]*$mps[4];
   $cnm   ++;
  };
  $count++; 
 };
 
 close(READER);
 
 if($cnm==1)
 {
  $m01  += $cm01;
  $m01d += $cm01d;
  $m23  += $cm23;
  $m23d += $cm23d;
  $j    += $cj;
  $jd   += $cjd;
  $j5   += $cj5;
  $j5d  += $cj5d;
  $nm   += $cnm;
  printf "\n\n\t FILE $file: ADDED TO STATISTICS, Q = %i!!!\n\n", $cnm;  
 };
 
};

printf "\n\n\t AVERAGES (over %i configs): \n\n", $nm;
printf "\t M01: %4.6lf +/- %4.6lf\n", $m01/$nm, sqrt(($m01d/$nm - ($m01*$m01)/($nm*$nm))/$nm);
printf "\t M23: %4.6lf +/- %4.6lf\n", $m23/$nm, sqrt(($m23d/$nm - ($m23*$m23)/($nm*$nm))/$nm);
printf "\t J:   %4.6lf +/- %4.6lf\n",   $j/$nm, sqrt(($jd/$nm - ($j*$j)/($nm*$nm))/$nm);
printf "\t J5:  %4.6lf +/- %4.6lf\n",  $j5/$nm, sqrt(($j5d/$nm - ($j5*$j5)/($nm*$nm))/$nm);

