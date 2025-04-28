#!/bin/bash

make

echo "Running ./q3_part1 "
./q3_part1 > q3_part1.out 2>&1

echo "Running ./q3_part2 -L 2" 
./q3_part2 -L 2 > q3_part2_lmax2.out 2>&1

echo "Running ./q3_part2 -L 8" 
./q3_part2 -L 8 > q3_part2_lmax8.out 2>&1

cd writeup
python3 plots.py
