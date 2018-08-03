# ./download.sh source destination
# means that you copy the file "source" from GRID to the local file "~/destination"

lcg-cp -D srmv2 -b -v --vo lattice.itep.ru srm://selattice.itep.ru:8446/srm/managerv2?SFN=/dpm/itep.ru/home/lattice.itep.ru/$LOGNAME/$1 file:$HOME/$2