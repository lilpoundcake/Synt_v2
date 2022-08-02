#!/bin/bash

apt-get update
apt-get install python3-dev
apt-get install gfortran
pip3 install --upgrade pip
pip3 install wheel
pip3 install numpy
pip3 install pandas
pip3 install primer3-py

git clone https://github.com/davidhoover/DNAWorks.git
cd DNAWorks/
make
cd ../Synt_v2

cp CHO.txt ../DNAWorks/
cp Sf9.txt ../DNAWorks/
cp bash_next.sh ../DNAWorks/
cp Logfile_generate.py ../DNAWorks/
cp Primer_Alg.py ../DNAWorks/
cp start.sh ../DNAWorks/

cd ../DNAWorks
bash start.sh
