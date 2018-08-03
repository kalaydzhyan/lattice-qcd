#!/usr/bin/perl
#$ -l h_rt=125:00:00

$bindir = sprintf("/afs/desy.de/user/t/tigrank/lattice/su3_over/bin");

while(true)
    {
	$cleanname = sprintf("%s/wd*H?",$bindir);
	@configs = glob("$cleanname");
	$nconfigs = scalar(@configs);
	if($nconfigs!=0)
	    {
		sleep($nconfigs*100);
	    	for($i=0; $i<$nconfigs; $i++)
	    	    {
	    		$config    = $configs[$i];
	    		system("rm -f $config");
	    		system("rm -f $config.dat");
	    	    }
	    }
	
	sleep(10);
    }