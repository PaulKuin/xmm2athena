#!/bin/csh
if ($1 == "HELLO") then
   echo This script from Paul is to add SUSS column where a match with Gaia/other 
   echo parameters: par1: dummy parameter   
endif
#
echo GO!
set IN1="/Users/data/catalogs/suss_gaia_epic/XMM-OM-SUSS5.0.fits#1"
set IN2="/Users/data/catalogs/gaia_xmm/id_merged_all.fit"
set OUT="/Users/data/catalogs/suss_gaia_epic/XMM-OM-SUSSA.fit"
#set SHORTOBSID=`echo $OBSID1 | bc`
echo "tmatch2 "
echo in1="/Users/data/catalogs/suss_gaia_epic/XMM-OM-SUSS5.0.fits#1" \
 in2=$IN2 out=$OUT
java -jar /Users/kuin/bin/topcat-full.jar -stilts tmatch2 \
#in1="/Users/data/catalogs/suss_gaia_epic/XMM-OM-SUSS5.0.fits#1" \
in1=$IN1 in2=$IN2 out=$OUT ifmt1=fits ifmt2=fits find=best1 join=all1 \
matcher=exact values1='IAUNAME' values2='IAUNAME' \
ocmd="addcol inGaiaeDR3 'length(IAUNAME_2) > 0';delcols 'IAUNAME_2 OBSID_2' "
