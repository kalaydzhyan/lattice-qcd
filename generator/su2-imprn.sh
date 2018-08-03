#PBS -N su28x4b2.95
#PBS -q long
#PBS -l cput=129:59:59,mem=900mb

#PBS -o /home/pools/2/lena/IMPROVED/L8x4/B2.950/out.out
#PBS -e /home/pools/2/lena/IMPROVED/L8x4/B2.950/err.err

cd /home/pools/2/lena/IMPROVED/L8x4/B2.950 

./su2-impr. 
exit 0
  
