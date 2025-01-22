#!/bin/csh 
#
#$HOME=/Users/kuin
set BIN=$HOME/bin
#
#  stilts_joined_mags3 - compute weighted value for three inputs to the magnitude
#
# if $0 < 10: not enough inputs
#  input and output file
set IN = $1
set OUT = $2
set NAME = $3
# input columns and their error in pairs (minumum 2)
set c4 = $4
set e5 = $5
set c6 = $6
set e7 = $7
set c8 = $8
set e9 = $9
echo "stilts_joined_mags3 in=$IN out=$OUT name=$NAME"
echo "columns $c4, $e5; $c6, $e7; $c8, $e9"
set EX = "($c4/($e5*$e5)+$c6/($e7*$e7)+$c8/($e9*$e9))/(1./($e5*$e5)+1./($e7*$e7)+1./($e9*$e9))"
echo "in shell : addcol $NAME $EX "
echo " "
java -jar $BIN/topcat-full.jar -stilts tpipe  in=$IN out=$OUT  ifmt=fits ofmt=fits cmd="addcol $NAME $EX "
  
