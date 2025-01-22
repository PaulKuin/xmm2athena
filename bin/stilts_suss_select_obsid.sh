#!/bin/csh
if ($2 == "") then
   echo This script from Paul is to extract the rows with a given OBSID from SUSS
   echo parameters: OBSID INFILE OUTFILE
endif
#alias stilts='java -jar /Users/kuin/bin/topcat-full.jar -stilts'
#echo $1
#echo $2
#echo $3
#
set OBSID1=$1
set IN=$2
set OUT=$3
set SHORTOBSID=`echo $OBSID1 | bc`
#echo java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe in=$IN out=$OUT cmd='select matches(OBSID,padWithZeros('${SHORTOBSID}',10))'
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe in=$IN ifmt=fits out=$OUT ofmt=fits \
cmd='select matches(OBSID,padWithZeros('${SHORTOBSID}',10))'

