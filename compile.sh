# CHANGE TO YOUR WORKING FOLDER!
cd /afs/desy.de/user/t/tigrank/lattice/su3_over

# if there exists gcc version 3.4 then we use it
# else we make a dummy soft-link to the current version of the gcc
# compiler in the current folder
rm -f ./gcc34
if [ -e /usr/bin/gcc34 ]
then
        ln -s -f /usr/bin/gcc34 ./gcc34
else
        ln -s -f /usr/bin/gcc ./gcc34
        gcc -v
fi

# here we're checking, whether there exists the old
# fortran compiler 
if [ -e /usr/bin/g77 ]
then
        cp var_old Makefile.var
else
        cp var_new Makefile.var
fi

# theese can be used in the case the gcc compiler can't
# find shared libraries like lapack in the ./lib folder

#export TEMPR=$LD_RUN_PATH
#export LD_RUN_PATH=$LD_RUN_PATH:$PWD/lib/
#export TEMPL=$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib/

cd 3dpart
make clean
make all
cd ..


# the first argument of this script = number of steps by x,y,z
# the second one = number of steps by t
# the third one = whether we need the single or double precision
# (1=single, 0=double)

if [ "$#" != "0"  ]
then
	make LAT_S=$1 LAT_T=$2 SINGLE=$3 clean
        make LAT_S=$1 LAT_T=$2 SINGLE=$3 all
        make LAT_S=$1 LAT_T=$2 SINGLE=$3 mcurr
else
	make clean
	make all
	make mcurr
fi

#export LD_RUN_PATH=$TEMPR
#export LD_LIBRARY_PATH=$TEMPL