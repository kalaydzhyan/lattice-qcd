#!/usr/bin/perl -w

$binary = "/home/pools/1/buividovich/overlap/bin/curr_1414s";

@Hs = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20);

foreach $H (@Hs)
{
### all H - beginning ###


$basedir = "/home/pools/1/buividovich/overlap/bin/b3.2810_s14_t14/H$H";

@files = glob("$basedir/ovd*.dat");

$nfiles = scalar(@files);

print "$nfiles files found.\n";

$maxnfiles = 50;
$nfiles = ($nfiles > $maxnfiles ? $maxnfiles : $nfiles);

print "Going to process $nfiles files...\n";

#Defining averages and dispersions

$maxq = 50;

@avX = (); @avX2 = ();
$nsamples = 0;

for($i=0; $i<$maxq; $i++)
{
 push(@avX, 0.0);
 push(@avX2, 0.0);
};

for($ifile=0; $ifile<$nfiles; $ifile++)
{
 print "Processing the file $files[$ifile]... \n";
 open(CURRENTS, "$binary -i $files[$ifile]|");
 $line = 0;
 while(<CURRENTS>)
 {
  print "$line: $_\n";
  @data = split;
  $ndata = scalar(@data);
  print "$ndata numbers obtained.\n";
  for($i=0; $i<$ndata; $i++)
  {
   $avX[$i]  += $data[$i];
   $avX2[$i] += $data[$i]*$data[$i];
  };
  $nsamples++;
  $line++;
 };
 close(CURRENTS);
};

for($i = 0; $i<$maxq; $i++)
{
 $avX[$i]  = $avX[$i]/$nsamples;
 $avX2[$i] = $avX2[$i]/$nsamples;
 $avX2[$i] = sqrt(($avX2[$i] - $avX[$i]*$avX[$i])/$nsamples);
};

@titles = ("ar5^2", "ar5*s23",
              "j0", 	         "j1",        	    "j2",	       "j3", 
            "j0^2",	       "j1^2",      	  "j2^2",      	     "j3^2",
        "r5^2j0^2",	   "r5^2j1^2",        "r5^2j2^2",        "r5^2j3^2", 
      "(d_0 r5)^2",      "(d_1 r5)^2",	    "(d_2 r5)^2",      "(d_3 r5)^2", 
             "s01",	 	"s02",             "s03",             "s12",   "s13",   "s23", 
           "s01^2",  	      "s02^2",           "s03^2",           "s12^2", "s13^2", "s23^2",
 "(d_0 r5)^2 j0^2", "(d_0 r5)^2 j1^2", "(d_0 r5)^2 j2^2", "(d_0 r5)^2 j3^2",
 "(d_1 r5)^2 j0^2", "(d_1 r5)^2 j1^2", "(d_1 r5)^2 j2^2", "(d_1 r5)^2 j3^2",
 "(d_2 r5)^2 j0^2", "(d_2 r5)^2 j1^2", "(d_2 r5)^2 j2^2", "(d_2 r5)^2 j3^2",
 "(d_3 r5)^2 j0^2", "(d_3 r5)^2 j1^2", "(d_3 r5)^2 j2^2", "(d_3 r5)^2 j3^2",
);
 
$ntitles = scalar(@titles) - 16; #Exclude dmar52jm2

print "\n\n\t STATISTICS: \t\n\n";

print "\n\t $nsamples samples \t\n";
 
for($i = 0; $i<$ntitles; $i++)
{
 printf "$titles[$i]:\t %6.6E \t\t +/- \t\t %6.6E, e = %3.2lf\n", $avX[$i], $avX2[$i], ($avX[$i]!=0)? abs($avX2[$i]/$avX[$i]) : 123.0;
};

print "\n\n\t CORRELATIONS: \t\n\n";

@cr52jm2 = ();

for($i=0; $i<4; $i++)
{
 $tmp = ($avX[10 + $i] - $avX[0]*$avX[6 + $i])/($nsamples*$avX2[0]*$avX2[6 + $i]);
 push(@cr52jm2, $tmp);
 printf "c(r5^2, j$i^2): %6.6E\n\n", $tmp;
 #$cr52jm2l1 = (($avX[10 + $i] - $avX2[10 + $i]) - ($avX[0] + $avX2[0])*($avX[6 + $i] + $avX2[6 + $i]))/($nsamples*$avX2[0]*$avX2[6 + $i]);
 #$cr52jm2l2 = (($avX[10 + $i] + $avX2[10 + $i]) - ($avX[0] - $avX2[0])*($avX[6 + $i] - $avX2[6 + $i]))/($nsamples*$avX2[0]*$avX2[6 + $i]);
 #$cr52jm2av = 0.5*($cr52jm2l1 + $cr52jm2l2);
 #$cr52jm2er = abs(0.5*($cr52jm2l1 - $cr52jm2l2));
 #printf "c(r5^2, j$i^2) = %6.6E +/- %6.6E\n\n", $cr52jm2av, $cr52jm2er; 
};

$cr5s23 = $avX[1]/sqrt($avX[0]*$avX[29]);

printf "c(r5, s23) = %6.6E\n\n", $cr5s23; 

#$cr5s23l1 = ($avX[1] - $avX2[1])/($avX[0]*$avX[25]);
#$cr5s23l2 = ($avX[1] + $avX2[1])/($avX[0]*$avX[25]);
#$cr5s23av = 0.5*($cr5s23l1 + $cr5s23l2);
#$cr5s23er = abs(0.5*($cr5s23l1 - $cr5s23l2));
#printf "c(r5, s23) = %6.6E +/- %6.6E\n\n", $cr5s23av, $cr5s23er; 

for($mu=0; $mu<4; $mu++)
{
 for($nu=0; $nu<4; $nu++)
 {
  $cdmar52jn2 = ($avX[30 + 4*$mu + $nu] - $avX[14 + $mu]*$avX[6 + $nu])/($nsamples*$avX2[6 + $nu]*$avX2[14 + $mu]);
  printf "%2.4E ", $cdmar52jn2;
 };
 print "\n";
}; 

open(DATA,">>cme.dat");

printf DATA "$H ";

printf DATA "%6.6E %6.6E ", $avX[0], $avX2[0]; # rho5^2


for($i=0; $i<4; $i++)
{
 printf DATA "%6.6E %6.6E ", $avX[6 + $i], $avX2[6 + $i]; # jm^2
};

for($i=0; $i<4; $i++)
{
 printf DATA "%6.6E %6.6E ", $avX[14 + $i], $avX2[14 + $i]; # (dm r5)^2
};

printf DATA "%6.6E %6.6E ", $avX[18], $avX2[18]; #s01
printf DATA "%6.6E %6.6E ", $avX[23], $avX2[23]; #s23
printf DATA "%6.6E %6.6E ", $avX[24], $avX2[24]; #s01^2
printf DATA "%6.6E %6.6E ", $avX[29], $avX2[29]; #s23^2

for($i=0; $i<4; $i++)
{
 printf DATA "%6.6E ", $cr52jm2[$i];
};

printf DATA "%6.6E ", $cr5s23;

printf DATA "\n";

close(DATA);


### All H - end ###
};