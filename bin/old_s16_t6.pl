#!/usr/bin/perl

$beta          = 3.3250;
$spacing       = 0.095;
$mass_phys     = 50.0; # IR Cutoff - in MeV
$max_nz_evals  = 10;    # UV cutoff
$LS            = 16;
$LT            = 6;
@Hs	       = (0, 1, 2, 3, 5, 10);

open(PARAMS, ">all.plt");

printf PARAMS "set term png\n";
printf PARAMS "set fit errorvariables\n";

$NH = scalar(@Hs);

for($j=1; $j<5; $j++)
{

 printf PARAMS "set print \"fitparam$j.txt\"\n";

 for($i=0; $i<$NH; $i++)
 {
  $H = $Hs[$i];

  printf PARAMS "f%i_%i(x)=a%i_%i + c%i_%i*(exp(-m%i_%i*x))\n", $H, $j, $H, $j, $H, $j, $H, $j;  
  printf PARAMS "fit f%i_%i(x) \"j_mu_nu_b%2.4f_H%i_s%i_t%i_m%i_ev%i.dat\" using 1:(-\$%i):(\$%i) via a%i_%i,c%i_%i,m%i_%i\n", $H, $j, $beta, $H, $LS, $LT, int($mass_phys), $max_nz_evals, 2*$j, 2*$j+1, $H, $j, $H, $j, $H, $j;
  printf PARAMS "print %i, a%i_%i, a%i_%i_err, m%i_%i, m%i_%i_err\n", $H, $H, $j, $H, $j, $H, $j, $H, $j;
  
 };

};

printf PARAMS "set key outside\n";
printf PARAMS "set key width +1\n";
printf PARAMS "set xlabel \"H\"\n";
printf PARAMS "set ylabel \"Mass parameter\"\n";

printf PARAMS "set output \"mass_s%i_t%i.png\"\n", $LS, $LT;

printf PARAMS "plot \"fitparam1.txt\" using 1:4:5 with yerrorbars title \"<j1 j1>\", ";
printf PARAMS "\"fitparam2.txt\" using 1:4:5 with yerrorbars title \"<j2 j2>\", ";
printf PARAMS "\"fitparam3.txt\" using 1:4:5 with yerrorbars title \"<j3 j3>\", ";
printf PARAMS "\"fitparam4.txt\" using 1:4:5 with yerrorbars title \"<j0 j0>\"\n";

printf PARAMS "set ylabel \"additive constant\"\n";

printf PARAMS "set output \"aconst_s%i_t%i.png\"\n", $LS, $LT;

printf PARAMS "plot \"fitparam1.txt\" using 1:2:3 with yerrorbars title \"<j1 j1>\", ";
printf PARAMS "\"fitparam2.txt\" using 1:2:3 with yerrorbars title \"<j2 j2>\", ";
printf PARAMS "\"fitparam3.txt\" using 1:2:3 with yerrorbars title \"<j3 j3>\", ";
printf PARAMS "\"fitparam4.txt\" using 1:2:3 with yerrorbars title \"<j0 j0>\"\n";

printf PARAMS "set logscale y\n";
printf PARAMS "set xlabel \"tau\"\n";

for($j=1; $j<5; $j++)
{

	if ($j==4) 
	{
		printf PARAMS "set output \"j0j0_s%i_t%i.png\"\n", $LS, $LT;
		printf PARAMS "set ylabel \"<j0 j0>\"\n";
	}
	else 
	{
		printf PARAMS "set output \"j%ij%i_s%i_t%i.png\"\n", $j, $j, $LS, $LT;
		printf PARAMS "set ylabel \"<j%i j%i>\"\n", $j, $j;
	};

	printf PARAMS "plot ";

	for($i=0; $i<$NH; $i++)
	{
	 $H = $Hs[$i];
	 printf PARAMS "f%i_%i(x) with lines title \"fit H=$H\", \"j_mu_nu_b%2.4f_H%i_s%i_t%i_m%i_ev%i.dat\" using 1:(-\$%i):%i with yerrorbars title \"H=$H\" ", $H, $j, $beta, $H, $LS, $LT, int($mass_phys), $max_nz_evals, 2*$j, 2*$j+1;
	 if ($i==($NH-1)) 
		{
		  printf PARAMS "\n";
		}
		 else
		 {
			printf PARAMS ", ";
		 };
	};

}

close(PARAMS);


system('gnuplot all.plt');
system('rm all.plt');