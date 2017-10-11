rm Run
rm -R output
mkdir output
cp src/Parameters.h output/Parameters.txt
#make clean
make
./Run
