cd ~/lattice/MC_GENERATOR/

if [ -f /usr/bin/gcc34 ]
then
	rm -f gcc34
	ln -s /usr/bin/gcc34 gcc34
else
	rm -f gcc34
        ln -s /usr/bin/gcc gcc34
        gcc -v
fi

#export TEMPR=$LD_RUN_PATH
#export LD_RUN_PATH=$LD_RUN_PATH:$PWD/lib/:$PWD/src/
#export TEMPL=$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/lib/:$PWD/src/

#echo $LD_LIBRARY_PATH

#cd 3dpart
#make clean
#make all
#cd ..

cd src
make clean
make all
cd ..

rm -f meanplaq
make meanplaq
rm -f alignment
make alignment
rm -f mc
make mc

#export LD_RUN_PATH=$TEMPR
#export LD_LIBRARY_PATH=$TEMPL

#ls -l /usr/lib

