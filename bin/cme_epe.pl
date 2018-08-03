#!/usr/bin/perl -w

$binary = "/home/pools/1/buividovich/overlap/bin/curr_1414s";

$basedir = "/home/pools/1/buividovich/overlap/bin/b3.2810_s14_t14/H10";

@files = glob("$basedir/ovd*Q+0.dat");

$nfiles = scalar(@files);

print "$nfiles files found.\n";

$maxnfiles = 20;
$nfiles = ($nfiles > $maxnfiles ? $maxnfiles : $nfiles);

print "Going to process $nfiles files...\n";

#Defining averages and dispersions

for($ifile=0; $ifile<$nfiles; $ifile++)
{
 print "Processing the file $files[$ifile]... \n";
 system("$binary -i $files[$ifile]");
};
