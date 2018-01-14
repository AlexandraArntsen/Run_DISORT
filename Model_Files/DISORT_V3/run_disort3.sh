#!/bin/sh
#
# This script compiles and executes DISORT
#

rm -f a.out
export OutputFile="disort_output.txt"
date > $OutputFile
cat DISOTEST3.f DISORT3.f BDREF.f DISOBRDF.f ERRPACK.f LINPAK.f LAPACK.f RDI1MACH.f > code.f
gfortran -O3 code.f 
chmod u+x ./a.out
time ./a.out
rm -f code.f
./a.out >> $OutputFile