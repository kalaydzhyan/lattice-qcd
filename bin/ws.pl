#!/usr/bin/perl -w

if ($#ARGV+1 < 7) {
    print "Usage: ws.pl TYPE R BETA SIZE_S SIZE_T NEV NCV_FACTOR_NOM [special]\n";
    print "       TYPE = {wd, ovd}\n";
    print "       special : if (TYPE == wd) { check ovd existance }\n";
    print "                 if (TYPE == ovd) { number of preconditioned modes }\n";
    die;
}

my $TYPE = $ARGV[0];
if (!($TYPE eq "wd" || "$TYPE" eq "ovd")) {
    die "TYPE should be either wd or ovd\n";
}
my $NEV_WD = 0;
my $CHECK_OVD = 0;
if ($TYPE eq "ovd") {
    if ($#ARGV+1 != 8) {
        die "For TYPE = ovd you should provide NEV_WD\n";
    }
    $NEV_WD = $ARGV[7];
} elsif ($TYPE eq "wd") {
    if ($#ARGV+1 == 8) {
        $CHECK_OVD = $ARGV[7];
    }
}

my $R      = $ARGV[1];
my $BETA   = $ARGV[2];
my $SIZE_S = $ARGV[3];
my $SIZE_T = $ARGV[4];
my $NEV    = $ARGV[5];
my $NCV_FACTOR_NOM = $ARGV[6];

my $DIR = "b${BETA}_s${SIZE_S}_t${SIZE_T}";
if (! -d $DIR ) {
    die "Directory $DIR does not exist";
}
if (! -x "./wd_${SIZE_S}${SIZE_T}s") {
    die "Wrong executable";
}

my $first = 1;
@lat_files = split( /\n/, `/bin/ls $DIR/00_b${BETA}_s${SIZE_S}_t${SIZE_T}_*.lat 2>/dev/null`);
@job_list = ("wd01","wd02","wd03","wd04","wd05","wd06","wd07","wd08","wd09","wd10","wd11","wd12","wd13","wd14","wd15","wd16",
"wd17","wd18","wd19","wd20","wd21","wd22","wd23","wd24","wd25","wd26","wd27","wd28","wd29","wd30","wd31","wd32","wd33","wd34","wd35","wd36","wd37","wd38");
$job_counter = 0;

print "Starting job scheduler\n";

for ( ; ; ) {
    sleep 4300 unless $first;
    if ($job_counter <= $#lat_files) {
        foreach (@job_list) {
            if ($job_counter > $#lat_files) { next; }
            @list = split(/\n/ , `qstat | grep smoroz | grep $_`);
            if (scalar @list == 0) {
                print "Attempt to start job '", $_ ,"' with data file ", $lat_files[$job_counter], "\n";
                if (run_next_job( $_, $lat_files[$job_counter]) == 0) {
                    $first = 0;
                }
                $job_counter++;
            }
        }
    } else {
        print "All jobs DONE, exiting\n";
        exit 0;
    }
}

###########################################################################
sub run_next_job {
    my $job = shift;
    my $data = shift;
    $data =~ s/$DIR\///;
    my @latfile = split(/_/, $data);
    my $number = $latfile[4];
    $number =~ s/\.lat$//;

    if ($TYPE eq "wd") {
        $wd_file = sprintf("wd_b%5.4f_s%i_t%i_%05i_sr_%2.1f_%i.dat",
                       ${BETA}, ${SIZE_S}, ${SIZE_T}, $number, $R, $NEV);
        if ($CHECK_OVD != 0) {
            $ovd_file = sprintf("ovd_%05i.r%2.1f.ev%i.dat", $number, $R, $CHECK_OVD);
            if ( -e "$DIR/$ovd_file" ) {
                return 1;
            }
        } else { 
      	    if ( -e "$DIR/$wd_file") {
                return 1;
            }
        }
    } elsif ($TYPE eq "ovd") {
        $wd_file = sprintf("wd_b%5.4f_s%i_t%i_%05i_sr_%2.1f_%i.dat",
                       ${BETA}, ${SIZE_S}, ${SIZE_T}, $number, $R, $NEV_WD);
        if ( !-e "$DIR/$wd_file") {
            return 0;
        }
        $ovd_file = sprintf("ovd_%05i.r%2.1f.ev%i.dat", $number, $R, $NEV);
        if ( -e "$DIR/$ovd_file" ) {
            return 1;
        }
    }
	



  
    open( FILE, "> run.sh" ) or die "Write open : $!";
    print FILE "#PBS -N $job\n";
    print FILE "#PBS -q long\n";
    print FILE "#PBS -l cput=199:59:59,mem=900mb\n\n";
    print FILE "cd /home/pools/2/lena/overlap/bin/$DIR\n";
    if ($TYPE eq "wd") {
        print FILE "../wd_${SIZE_S}${SIZE_T}s -i lwd_$number.info -o $wd_file -h -s -r $R -S $NEV -L 2 -p -m 0.5 -d 15 -a -F $NCV_FACTOR_NOM\n";
    } elsif ($TYPE eq "ovd") {
        print FILE "../ovd_${SIZE_S}${SIZE_T}s -i lwd_$number.info -O $wd_file -o $ovd_file -r $R -S $NEV -t 0 -a\n";
    }
    print FILE "exit 0\n";
    close FILE;
 
    die "Qsub got fatal error : $!" unless system("run.sh", "run.sh") == 0;
    sleep 2;
    return 0;
}

# vim: ts=4 : sw=4 : expandtab :
