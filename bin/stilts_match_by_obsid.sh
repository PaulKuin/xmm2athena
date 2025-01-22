#!/bin/csh
if ($2 == "") then
   echo This script from Paul is to match SUSS rows to Gaia for fixed OBSID
   echo parameters: INFILE1: SUSS_obsid INFILE2: Gaia_obsid OUTFILE
   echo example:   
   echo   stilts tmatch2 in1=suss_0000110101.fits in2=gaiaedr3_for_0000110101_a.vot \
   echo   out=suss5.0_gaiaedr3_0000110101.fits  matcher=sky values1="RA DEC" 
   echo   values2="ra dec" params=1.0 find=best1
endif
set IN1 = $1
set IN2 = $2
set OUT = $3
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmatch2 in1=$IN1 in2=$IN2 out=$OUT matcher=sky values1="RA DEC" values2="ra dec" params=1.0 find=best1
