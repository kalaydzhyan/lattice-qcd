#!/usr/bin/perl

 $bindir = sprintf("./");
 
 
 $cleanname = sprintf("%s/ttt*",$bindir);
 @cores = glob("$cleanname");
 $ncores = scalar(@cores);
 if($ncores!=0)
 {
 system("ls -l $bindir/*.e* > ls.dat");
 system("mail -s \'Shit happens...\' kalaydzhyan\@gmail.com < ls.dat");
 
 }