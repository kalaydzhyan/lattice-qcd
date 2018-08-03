#!/usr/bin/perl

 @m = glob("wait*");
 $n = scalar(@m);
 
 if($#ARGV==1)
 {
   $tname = sprintf("run%03is_H%i.sh", $ARGV[0],$ARGV[1]);
   
   if($n==0) 
     {
        system("qsub -l arch=x86 -l h_vmem=2G -cwd b8.3000_s14_t14/$tname");
#        system("qsub -l -cwd b8.4500_s16_t6/$tname");
     }
   else
     {
        print("please wait and try again...\n");
     }
 }