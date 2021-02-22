#!/usr/bin/env bash

cd $toolsDir
mkdir skesa
cd skesa
git clone https://github.com/ncbi/SKESA.git
cd SKESA
conda install -c conda-forge/label/gcc7 boost
make -f Makefile.nongs
