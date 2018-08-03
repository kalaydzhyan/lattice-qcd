#!/usr/bin/perl

$BETABEGIN  =  6.0;
$BETAEND    =  9.9;

for($i=$BETABEGIN; $i<=$BETAEND; $i+=0.1)
{

 $scriptname = sprintf("beta_%2.3f.sh", $i);

 system("rm -f $scriptname");

 open(PARAMS, ">$scriptname");
# printf PARAMS "#PBS -N beta_%2.3f\n#PBS -q long\n#PBS -l cput=199:59:59,mem=900mb\n", $i;
 printf PARAMS "cd ~/lattice/MC_GENERATOR/\n";
 printf PARAMS "./alignment %f\n", $i;
 close(PARAMS);
 
 system("qsub -cwd $scriptname");
 print "batch for beta=$i has been added\n";

};