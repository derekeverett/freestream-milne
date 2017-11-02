export OMP_NUM_THREADS=4 #set this to the number of threads on CPU
rm Run
rm -R output
mkdir output
cp src/Parameters.h output/Parameters.txt
#make clean
make
./Run
