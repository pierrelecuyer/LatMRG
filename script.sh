#!/bin/sh
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
