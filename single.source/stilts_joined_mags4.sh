#!/bin/csh 
#
#  stilts_joined_mags4 - compute weighted value for four inputs to the magnitude
#
#$HOME=/Users/kuin
set BIN=$HOME/bin
#
# if $0 < 12: not enough inputs
echo "stilts_joined_mags4 in=$1 out=$2 name=$3"
echo "columns $4, $5; $6, $7; $8,$9; $10, $11"
#  input and output file
$IN = $1
$OUT = $2
$NAME = $3
# input columns and their error in pairs (minumum 2)
# $col1 = $4
# $err1 = $5
# $col2 = $6
# $err2 = $7
# $col3 = $8
# $err3 = $9
# $col4 = $10
# $err4 = $11
java -jar $BIN/topcat-full.jar -stilts tpipe \
  in=$IN \
  out=$OUT  \
  ifmt='fits' ofmt='fits' \ 
  cmd="addcol $NAME '($4/($5*$5) + $6/($7*$7) + $8/($9*$9)+ $10/($11*$11))/(1./($5*$5)+1./($7*$7)+1./($9*$9)+1./($11*$11))' "
