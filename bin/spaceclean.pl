#!/usr/bin/perl
#$ -l h_rt=125:00:00

$bindir = sprintf("/afs/desy.de/user/t/tigrank/lattice/su3_over/bin");

while(true)
    {
	$cleanname = sprintf("%s/core*",$bindir);
	@cores = glob("$cleanname");
	$ncores = scalar(@cores);
	 if($ncores!=0)
	    {
	       system("ls -l $bindir/*.e* > ls.dat");
	       system("mail -s \'Shit happens...\' kalaydzhyan\@gmail.com < ls.dat");
	       system("rm $bindir/ls.dat");
    	    }

	system("rm -f $bindir/core*");
	
#	$cleanname = sprintf("%s/wd*.dat",$bindir);
#	@configs = glob("$cleanname");
#	$nconfigs = scalar(@configs);
#	if($nconfigs!=0)
#	    {
#		sleep($nconfigs*600);
#	    	for($i=0; $i<$nconfigs; $i++)
#	    	    {
#	    		$config    = $configs[$i];
#	    		system("rm -f $config");
#	    	    }
#	    }
	
	sleep(10);
    }