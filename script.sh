#!/bin/sh
# This is a working script. This script uses sleep because yafu (or temporary files
# for factorization) cannot be instanced multiple times at once
./bin/SeekRe examples/inputs/seekre/pow2/pow2_31_31_3.dat > 31_31_3.res &
sleep 600
./bin/SeekRe examples/inputs/seekre/pow2/pow2_31_41_4.dat > 31_41_4.res &
sleep 600
./bin/SeekRe examples/inputs/seekre/pow2/pow2_31_51_5.dat > 31_51_5.res &
sleep 600
./bin/SeekRe examples/inputs/seekre/pow2/pow2_61_31_3.dat > 61_31_3.res &
sleep 600
./bin/SeekRe examples/inputs/seekre/pow2/pow2_61_41_4.dat > 61_41_4.res &
sleep 600
./bin/SeekRe examples/inputs/seekre/pow2/pow2_61_51_5.dat > 61_51_5.res &
