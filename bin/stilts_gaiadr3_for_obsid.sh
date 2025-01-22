#!/bin/csh
#  {ra_} {dec_} {radius} {gaia-infile} {gtemp}
if ($2 == "") then
   echo This script from Paul is to extract the cone with radius from high PM gaia dr3
   echo parameters: obsid ra dec radius gaia-infile gaiatempfile
endif
#
set RA=$1
set DEC=$2
set RADIUS=$3
set IN=$4
set OUT=$5
echo pointing $RA, $DEC
echo radius  $RADIUS 
echo input file $IN
echo output file $OUT
# select cone from input catalogue
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe in=$IN ifmt=fits out=$OUT ofmt=fits \
cmd="select skyDistance(ra,dec,$RA,$DEC)<$RADIUS"
#

