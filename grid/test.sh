tar -xvf toverlap1.tar
cd toverlap
./compile.sh 4 4 1
cd bin
./run.pl
echo 
echo 
cat out.dat
