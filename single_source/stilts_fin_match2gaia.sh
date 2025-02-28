#!/bin/csh 
#
echo "finalise.match2Gaia.sh: rematch all"
set IN1 = matchout_all.fits  
set IN2 = nomatch_all.fits
set OUT = sussxgaiadr3_ep2000.fits
#
echo first file
echo
java -jar ~/bin/topcat-full.jar -stilts cdsskymatch cdstable=I/355/gaiadr3 \
 find=each in=$IN1  ifmt=fits  ra=ra2000Ep dec=dec2000Ep radius=3 out=gaia_pm_match.fits \
 ocmd='delcols "ra_error dec_error parallax_over_error healpix raObs decObs obsEp 2000Ep" '\
 ocmd='addcol mismatchedPM "((pm_in - pm_cds > 5)?1:0)" '
echo 
echo second file  
echo
java -jar ~/bin/topcat-full.jar -stilts cdsskymatch cdstable=I/355/gaiadr3 find=each in=$IN2 ifmt=fits \
 ra=ra2000Ep dec=dec2000Ep radius=3 out=gaia_no_pm_match.fits \
 ocmd='addcol mismatchedPM "0" ' ocmd='delcols "2000Ep"'
echo
echo make file 2 match preparing
java -jar ~/bin/topcat-full.jar -stilts tpipe in=gaia_pm_match.fits ifmt=fits  out=tmp_gaia_pm_match.fits ofmt=fits cmd='delcols "designation source_id ra2016gaia dec2016gaia parallax pm_in pmra_in pmra_error pmdec pmdec_error Separation"' 
echo
echo tcat them
java -jar ~/bin/topcat-full.jar -stilts tcat in='tmp_gaia_pm_match.fits gaia_no_pm_match.fits' ifmt=fits out=$OUT ofmt=fits
echo
