#!/usr/bin/perl -w

$binary = "/home/pools/1/buividovich/overlap/bin/tpc_1616s";

$basedir    = "/home/pools/2/lena/overlap_16x16_ext/bin/b3.3555_s16_t16";
$targetdir  = "/home/pools/1/buividovich/overlap/bin/b3.3555_s16_t16";

@dirs = glob("$basedir/field_H=*_E=0/");
$ndirs = scalar(@dirs);
print "$ndirs directories found ...\n";

foreach $dir (@dirs)
{
 # Extracting the magnetic field and making the new directory in the target folder
 @pps = split('H=',$dir);
 chop($pps[1]);
 system("mkdir $targetdir/H$pps[1]");
 # Processing and copying
 print "\t Processing and copying the directory $dir ... \n";
 @files = glob("$dir/ovd*.dat");
 $nfiles = scalar(@files);
 print "$nfiles files found...\n";
 $filecount = 1;
 foreach $file (@files)
 {
  print "\t \t Processing and copying the file $file ...\n";
  open(TPC, "$binary -i $file|");
  while(<TPC>)
  {
   $Q = $_;
  };
  close(TPC);
  chop($Q);
  print "\t\t Topological charge: $Q\n";
  $newfname = sprintf("$targetdir/H$pps[1]/ovd_%05i_ev30_Q%+i.dat", $filecount, $Q);
  print "Copying to: $newfname\n";
  system("cp -s $file $newfname");
  $filecount ++;
 };
};
