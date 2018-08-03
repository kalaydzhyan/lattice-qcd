# ./upload.sh source destination
# means that you copy the local file "source" the GRID file "~/destination"

lcg-cp -D srmv2 -b -v --vo lattice.itep.ru file:$HOME/$1 srm://selattice.itep.ru:8446/srm/managerv2?SFN=/dpm/itep.ru/home/lattice.itep.ru/$LOGNAME/$2