#!/bin/bash

./overlapBed.pl -i1 Intron.bed -i2 TE.bed -m 50 -j -o output # overlap Intron.bed and TE.bed, print out the pairs of overlapping entries from both bed files that have minimum percent overlap of 50% 


###############################################
##
## Your code should run with, BUT NOT LIMITED TO, the command(s) above.
## Please read the exercise instructions for more details on how your code should work.
##
###############################################