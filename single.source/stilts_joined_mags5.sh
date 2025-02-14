#!/bin/csh 
#
#  stilts_joined_mags5 - compute weighted value for five inputs to the magnitude
#
#$HOME=/Users/kuin
set BIN=$HOME/bin
#
#
# if $0 < 14: not enough inputs
echo "stilts_joined_mags5 in=$1 out=$2 name=$3"
echo "columns $4, $5; $6, $7; $8,$9; $10,$11; $11,$12"
#  input and output file
set IN = $1
set OUT = $2
set NAME = $3
# input columns and their error in pairs (minumum 2)
# $col1 = $4
# $err1 = $5
# $col2 = $6
# $err2 = $7
# $col3 = $8
# $err3 = $9
# $col4 = $10
# $err4 = $11
# $col5 = $12
# $err5 = $13
java -jar $BIN/topcat-full.jar -stilts tpipe \
  in=$IN \
  out=$OUT  \
  ifmt='fits' ofmt='fits' \ 
  cmd="addcol $NAME '($4/($5*$5) + $6/($7*$7) + $8/($9*$9)+ $10/($11*$11) + $12/($13*$13))/(1./($5*$5)+1./($7*$7)+1./($9*$9)+1./($11*$11)+1./($13*$13))' "
  
