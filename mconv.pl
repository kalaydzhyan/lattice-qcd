#!/usr/bin/perl

@Hs = (1, 2, 3, 5, 10);

foreach $H (@Hs)
{
 open(DATA, "</home/pools/2/lena/overlap_16x6_ext/bin/b3.1600_s16_t6/field_H=$H/data.txt");
 open(MDATA, ">/home/pools/1/buividovich/overlap/mdata_16x6_3.1600_H$H.dat");
 while(<DATA>)
 {
  $_ =~ s/[(),]/ /g;
  @fields = split(/\s+/);
  $nfields= scalar(@fields);
  printf "\nSplitted into %i parts\n", scalar(@fields);
  print       "$fields[1] $fields[2] 0 $fields[3] $fields[4]\n";
  print MDATA "$fields[1] $fields[2] 0 $fields[3] $fields[4]\n";
 };
 close(MDATA);
 close(DATA);
};