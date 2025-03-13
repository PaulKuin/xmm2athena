#!/bin/csh 
#
#$HOME=/Users/kuin
set BIN=$HOME/bin
#
#  stilts_joined_mags2 - compute weighted value for two inputs to the magnitude
#
# if $0 < 8: not enough inputs
echo "stilts_joined_mags2 in=$1 out=$2 name=$3"
echo "columns $4, $5; $6, $7"
#  input and output file
set IN = $1
set OUT = $2
set NAME = $3
# input columns and their error in pairs (minumum 2)
# $col1 = $4
# $err1 = $5
# $col2 = $6
# $err2 = $7
#
java -jar $BIN/topcat-full.jar -stilts tpipe in=$1  out=$2  \
  ifmt='fits' ofmt='fits' \ 
  cmd="addcol $3  ($4/($5*$5) + $6/($7*$7))/(1./($5*$5)+1./($7*$7))"
