#!/usr/bin/perl -w

$mssc = "/home/pools/1/buividovich/overlap/bin/mssc_166s";
#$csp  = 0.226044; # beta = 3.2810
$csp   = 0.277296; # beta = 3.1600
#$csp  = 0.208859; #beta = 3.3250

die "Correct usage: magnetization.pl indir\n" unless (scalar(@ARGV)==1); 

$indir = $ARGV[0];

@flist = glob "$indir/ovd*.dat";
$nfiles = scalar(@flist);

$maxnfiles = ($nfiles > 20)? 20 : $nfiles;

print "\n\n\t $nfiles files found, processing the first $maxnfiles only.\n\n";

#determining the topological charges of the configurations

$nbins    = 10;
$binsize  = 50;
@bins_m01   = ();
@bins_m01_2 = ();
@bins_num = ();

for($i=0; $i<$nbins; $i++)
{
 push(@bins_m01,   0.0);
 push(@bins_m01_2, 0.0);
 push(@bins_num,     0);
};  

for($ifile=0; $ifile<$maxnfiles; $ifile++)
{
 $file = $flist[$ifile];
 print "Analyzing the file $file ... \n";
 
 open(READER, "$mssc -i $file -c $csp|");
 $count  = 0;
 while(<READER>)
 {
  @mps = split;
  die "Wrong output of mssc!!!\n" unless (scalar(@mps)==2);
  $ibin = int($mps[0]/$binsize);
  if($mps[0] > 1)
  {
   $bins_m01[$ibin]   += $mps[1];
   $bins_m01_2[$ibin] += $mps[1]*$mps[1];
   $bins_num[$ibin] ++;
  };    
  $count++; 
 };
 close(READER);
};

for($i=0; $i<$nbins; $i++)
{
 if($bins_num[$i]!=0)
 {
  $bins_m01[$i]   = $bins_m01[$i]  /$bins_num[$i];
  $bins_m01_2[$i] = $bins_m01_2[$i]/$bins_num[$i];
  $bins_m01_2[$i] = 1/sqrt($bins_num[$i])*sqrt($bins_m01_2[$i] - $bins_m01[$i]*$bins_m01[$i]);
 }; 
 printf "%i: %i data points, m01 =  %2.4f +/- %2.4f\n", $i, $bins_num[$i], $bins_m01[$i], $bins_m01_2[$i];
};
