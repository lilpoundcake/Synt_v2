#!/bin/bash


#/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
#echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zprofile
#eval "$(/opt/homebrew/bin/brew shellenv)"
git clone https://github.com/Homebrew/brew homebrew
eval "$(homebrew/bin/brew shellenv)"
brew update --force --quiet
chmod -R go-w "$(brew --prefix)/share/zsh"

brew install python@3.10
brew install gcc
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
