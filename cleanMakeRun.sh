export OMP_NUM_THREADS=40 #set this to the number of threads on CPU
rm Run
rm -R output
mkdir output
cp parameters.dat output/parameters.dat
#make clean
make
./Run
