#!/bin/bash

./dnaworks р_1.inp
python3 Primer_Alg.py
mkdir р
cp р_1.inp р
cp р_1.txt р
cp р_hairpin.csv р
cp р_high_temp.csv р
cp р_SG_primers.fasta р
cp р_sequence.fasta р
cd
cp -r ./DNAWorks/р/ .
rm -rf ./DNAWorks ./For_DNAWorks
tar -zcvf р.tar.gz ./р/
