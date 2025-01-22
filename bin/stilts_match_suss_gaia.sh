#!/bin/csh
if ($1 == "") then
   echo This script from Paul is to match SUSS rows to Gaia for fixed OBSID
   echo parameters: INFILE1: SUSS_obsid INFILE2: Gaia_obsid OUTFILE
   echo example:   
   echo   stilts tmatch2 in1=suss_0000110101.fits in2=gaiaedr3_for_0000110101_a.vot \
   echo   out=suss5.0_gaiaedr3_0000110101.fits  matcher=sky values1="RA DEC" 
   echo   values2="ra dec" params=1.0 find=best1
   echo This script from Paul is to extract the rows with a given OBSID from SUSS
   echo parameters: OBSID  OUTFILE
endif
#
set OBSID1=$1
#set IN=\"/Users/data/catalogs/suss_gaia_epic/XMM-OM-SUSS5.0.fits\#1\"
set TMP=/Users/data/catalogs/gaia_xmm/matched/temp.fits
set SHORTOBSID=`echo $OBSID1 | bc`
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe \
in="/Users/data/catalogs/suss_gaia_epic/XMM-OM-SUSS5.0.fits#1" \
out=$TMP cmd='select matches(OBSID,padWithZeros('${SHORTOBSID}',10))'
#
set IN1 = $TMP
set IN2 = "/Users/data/catalogs/gaia_xmm/adir/gaiaedr3_for_"$OBSID1"_a.vot"
set OUT = "/Users/data/catalogs/gaia_xmm/matched/suss5.0_gaiaedr3_"$OBSID1".fits"
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmatch2 in1=$IN1 in2=$IN2 \
out=$OUT matcher=sky values1="RA DEC" values2="ra dec" params=1.0 find=best1 join=all1
