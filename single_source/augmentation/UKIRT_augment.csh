alias topcat 'java -jar /Users/kuin/bin/topcat-full.jar '
# revision npmk 2025-03-09
echo 
echo Add petrosian photometry to CDS-upload-xmatch UKIDSS datasets  
date
echo
echo now DXS survey
topcat -stilts coneskymatch in=SUSS_DXS.fits ifmt=fits out=out.fits icmd='select (ra>0.0)' find=each ra=ra2000Ep dec=dec2000Ep servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/dxs9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra dec pJmag e_pJmag pKmag e_pKmag";'
#
topcat -stilts tskymatch2 in1=SUSS_DXS.fits ifmt1=fits in2=out.fits out=out2.fits ra1='ra2000Ep' dec1='dec2000Ep' ra2='ra' dec2='dec' error=0.5 find=best1
topcat -stilts tpipe in=out2.fits out=SUSS_DXS_pet.fits cmd='delcols "ra_2 dec_2 Separation"'
date
#
echo
echo now GCS survey takes 1 hr to run
echo
topcat -stilts coneskymatch in=SUSS_GCS.fits ifmt=fits out=out3.fits icmd='select (ra>0.0)' find=each ra=ra2000Ep dec=dec2000Ep servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/gcs9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='keepcols "ra dec pZmag e_pZmag pYmag e_pYmag pJmag e_pJmag pHmag e_pHmag pKmag1 e_pKmag1";'
#
topcat -stilts tskymatch2 in1=SUSS_GCS.fits ifmt1=fits in2=out3.fits out=out4.fits ra1='ra2000Ep' dec1='dec2000Ep' ra2='ra' dec2='dec' error=0.5 find=best1
topcat -stilts tpipe in=out4.fits out=SUSS_GCS_pet.fits cmd='delcols "ra_2 dec_2 Separation"'
date
echo
echo
echo  the LAS processing takes 4-5 hrs to run
echo
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out1.fits icmd='select (ra>0.0)' find=each ra=ra2000Ep dec=dec2000Ep servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='keepcols "ra dec pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";'
#
date
#topcat -stilts tcat in="out1.fits out2.fits out3.fits out4.fits out5.fits out6.fits out7.fits out8.fits out10.fits out11.fits" out=out.fits 
topcat -stilts tskymatch2 in1=SUSS_LAS.fits ifmt1=fits in2=out1.fits out=out.fits ra1='ra2000Ep' dec1='dec2000Ep' ra2='ra' dec2='dec' error=0.5 find=best1
topcat -stilts tpipe in=out.fits out=SUSS_LAS_pet.fits cmd='delcols "ra_2 dec_2 GroupID GroupSize Separation"'
date
echo
echo Forget about GPS survey as there are no petrosian magnitudes in GPS.
echo
echo
echo Now connect all the UKIDSS surveys together
echo
topcat -stilts tmatch2 in1=SUSS_GCS_pet.fits ifmt1=fits in2=SUSS_DXS_pet.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd1=' colmeta -name GCS_zmag "zAperMag3"; colmeta -name GCS_e_zmag "zAperMag3Err"; colmeta -name GCS_ymag "yAperMag3"; colmeta -name GCS_e_ymag "yAperMag3Err"; colmeta -name GCS_jmag "jAperMag3"; colmeta -name GCS_e_jmag "jAperMag3Err"; colmeta -name GCS_hmag "hAperMag3"; colmeta -name GCS_e_hmag "hAperMag3Err"; colmeta -name GCS_kmag "k_1AperMag3"; colmeta -name GCS_e_kmag "k_1AperMag3Err"; colmeta -name GCS_zPmag "pZmag"; colmeta -name GCS_e_zPmag "e_pZmag"; colmeta -name GCS_yPmag "pYmag"; colmeta -name GCS_e_yPmag "e_pYmag"; colmeta -name GCS_jPmag "pJmag"; colmeta -name GCS_e_jPmag "e_pJmag"; colmeta -name GCS_hPmag "pHmag"; colmeta -name GCS_e_hPmag "e_pHmag"; colmeta -name GCS_kPmag "pKmag1"; colmeta -name GCS_e_kPmag "e_pKmag1"; delcols "JName ra_1 dec_1 sourceID mode epoch mergedClass angDist" ' icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; colmeta -name DXS_jmag "jAperMag3"; colmeta -name DXS_e_jmag "jAperMag3Err"; colmeta -name DXS_kmag "kAperMag3"; colmeta -name DXS_e_kmag "kAperMag3Err"; colmeta -name DXS_jPmag "pJmag"; colmeta -name DXS_e_jPmag "e_pJmag"; colmeta -name DXS_kPmag "pKmag"; colmeta -name DXS_e_kPmag "e_pKmag"; keepcols "SRCNUM2 OBSID2 DXS_jmag DXS_e_jmag DXS_kmag DXS_e_kmag DXS_jPmag DXS_e_jPmag DXS_kPmag DXS_e_kPmag"' ocmd='delcols "SRCNUM2 OBSID2 Separation"'
cp temp.fits temp_in.fits
echo 
echo add GPS survey
date
echo 
topcat -stilts tmatch2 in1=temp_in.fits ifmt1=fits in2=SUSS_GPS.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; colmeta -name GPS_jmag "jAperMag3"; colmeta -name GPS_e_jmag "jAperMag3Err"; colmeta -name GPS_hmag "hAperMag3"; colmeta -name GPS_e_hmag "hAperMag3Err"; colmeta -name GPS_kmag "k_1AperMag3"; colmeta -name GPS_e_kmag "k_1AperMag3Err"; colmeta -name GP2_kmag "k_2AperMag3"; colmeta -name GP2_e_kmag "k_2AperMag3Err"; keepcols "SRCNUM2 OBSID2 GPS_jmag GPS_e_jmag GPS_hmag GPS_e_hmag GPS_kmag GPS_e_kmag GP2_kmag GP2_e_kmag"' ocmd='delcols "SRCNUM2 OBSID2 GroupID GroupSize Separation"'
cp temp.fits temp_in.fits
echo 
echo add LAS survey
date
echo 
topcat -stilts tmatch2 in1=temp_in.fits ifmt1=fits in2=SUSS_LAS_pet.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; colmeta -name LAS_ymag "yAperMag3"; colmeta -name LAS_e_ymag "yAperMag3Err"; colmeta -name LAS_jmag "j_1AperMag3"; colmeta -name LAS_e_jmag "j_1AperMag3Err"; colmeta -name LAS_hmag "hAperMag3"; colmeta -name LAS_e_hmag "hAperMag3Err"; colmeta -name LAS_kmag "kAperMag3"; colmeta -name LAS_e_kmag "kAperMag3Err"; colmeta -name LAS_yPmag "pYmag"; colmeta -name LAS_e_yPmag "e_pYmag"; colmeta -name LAS_jPmag "pJmag1"; colmeta -name LAS_e_jPmag "e_pJmag1"; colmeta -name LAS_hPmag "pHmag"; colmeta -name LAS_e_hPmag "e_pHmag"; colmeta -name LAS_kPmag "pKmag"; colmeta -name LAS_e_kPmag "e_pKmag"; keepcols "SRCNUM2 OBSID2 LAS_ymag LAS_e_ymag LAS_jmag LAS_e_jmag LAS_hmag LAS_e_hmag LAS_kmag LAS_e_kmag LAS_yPmag LAS_e_yPmag LAS_jPmag LAS_e_jPmag LAS_hPmag LAS_e_hPmag LAS_kPmag LAS_e_kPmag"' ocmd='delcols "SRCNUM2 OBSID2 Separation"'
cp temp.fits temp_in.fits
echo
echo Now onto magnitude prioritization
echo 
echo Prioritize K magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits cmd='addcol -units mag UK_kmag "( (GPS_kmag > 0.0) ? GPS_kmag : GP2_kmag )";addcol -units mag UK_e_kmag "( (GPS_kmag > 0.0) ? GPS_e_kmag : GP2_e_kmag )";delcols "GPS_kmag GPS_e_kmag GP2_kmag GP2_e_kmag";colmeta -name GPS_kmag "UK_kmag";colmeta -name GPS_e_kmag "UK_e_kmag";addcol -units mag UK_kmag "( (GPS_kmag > 0.0) ? GPS_kmag : LAS_kmag )";addcol -units mag UK_e_kmag "( (GPS_kmag > 0.0) ? GPS_e_kmag : LAS_e_kmag )";delcols "GPS_kmag GPS_e_kmag LAS_kmag LAS_e_kmag";colmeta -name GPS_kmag "UK_kmag";colmeta -name GPS_e_kmag "UK_e_kmag";addcol -units mag UK_kmag "( (DXS_kmag > 0.0) ? DXS_kmag : GPS_kmag )";addcol -units mag UK_e_kmag "( (DXS_kmag > 0.0) ? DXS_e_kmag : GPS_e_kmag )";delcols "DXS_kmag DXS_e_kmag GPS_kmag GPS_e_kmag";colmeta -name DXS_kmag "UK_kmag";colmeta -name DXS_e_kmag "UK_e_kmag";addcol -units mag UK_kmag "( (GCS_kmag > 0.0) ? GCS_kmag : DXS_kmag )";addcol -units mag UK_e_kmag "( (GCS_kmag > 0.0) ? GCS_e_kmag : DXS_e_kmag )";delcols "GCS_kmag GCS_e_kmag DXS_kmag DXS_e_kmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize H magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits cmd='addcol -units mag UK_hmag "( (GPS_hmag > 0.0) ? GPS_hmag : LAS_hmag )";addcol -units mag UK_e_hmag "( (GPS_hmag > 0.0) ? GPS_e_hmag : LAS_e_hmag )";delcols "GPS_hmag GPS_e_hmag LAS_hmag LAS_e_hmag";colmeta -name GPS_hmag "UK_hmag";colmeta -name GPS_e_hmag "UK_e_hmag";addcol -units mag UK_hmag "( (GCS_hmag > 0.0) ? GCS_hmag : GPS_hmag )";addcol -units mag UK_e_hmag "( (GCS_hmag > 0.0) ? GCS_e_hmag : GPS_e_hmag )";delcols "GCS_hmag GCS_e_hmag GPS_hmag GPS_e_hmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize J magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits cmd='addcol -units mag UK_jmag "( (GPS_jmag > 0.0) ? GPS_jmag : LAS_jmag )";addcol -units mag UK_e_jmag "( (GPS_jmag > 0.0) ? GPS_e_jmag : LAS_e_jmag )";delcols "GPS_jmag GPS_e_jmag LAS_jmag LAS_e_jmag";colmeta -name GPS_jmag "UK_jmag";colmeta -name GPS_e_jmag "UK_e_jmag";addcol -units mag UK_jmag "( (DXS_jmag > 0.0) ? DXS_jmag : GPS_jmag )";addcol -units mag UK_e_jmag "( (DXS_jmag > 0.0) ? DXS_e_jmag : GPS_e_jmag )";delcols "DXS_jmag DXS_e_jmag GPS_jmag GPS_e_jmag";colmeta -name DXS_jmag "UK_jmag";colmeta -name DXS_e_jmag "UK_e_jmag";addcol -units mag UK_jmag "( (GCS_jmag > 0.0) ? GCS_jmag : DXS_jmag )";addcol -units mag UK_e_jmag "( (GCS_jmag > 0.0) ? GCS_e_jmag : DXS_e_jmag )";delcols "GCS_jmag GCS_e_jmag DXS_jmag DXS_e_jmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize Y magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits cmd='addcol -units mag UK_ymag "( (GCS_ymag > 0.0) ? GCS_ymag : LAS_ymag )";addcol -units mag UK_e_ymag "( (GCS_ymag > 0.0) ? GCS_e_ymag : LAS_e_ymag )";delcols "GCS_ymag GCS_e_ymag LAS_ymag LAS_e_ymag"'
cp temp.fits temp_in.fits
echo 
echo Rename Z magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits cmd='colmeta -name UK_zmag "GCS_zmag";colmeta -name UK_e_zmag "GCS_e_zmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize Petrosian K magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits cmd='addcol -units mag UK_kPmag "( (DXS_kPmag > 0.0) ? DXS_kPmag : LAS_kPmag )";addcol -units mag UK_e_kPmag "( (DXS_kPmag > 0.0) ? DXS_e_kPmag : LAS_e_kPmag )";delcols "DXS_kPmag DXS_e_kPmag LAS_kPmag LAS_e_kPmag";colmeta -name DXS_kPmag "UK_kPmag";colmeta -name DXS_e_kPmag "UK_e_kPmag";addcol -units mag UK_kPmag "( (GCS_kPmag > 0.0) ? GCS_kPmag : DXS_kPmag )";addcol -units mag UK_e_kPmag "( (GCS_kPmag > 0.0) ? GCS_e_kPmag : DXS_e_kPmag )";delcols "GCS_kPmag GCS_e_kPmag DXS_kPmag DXS_e_kPmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize Petrosian H magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits cmd='addcol -units mag UK_hPmag "( (GCS_hPmag > 0.0) ? GCS_hPmag : LAS_hPmag )";addcol -units mag UK_e_hPmag "( (GCS_hPmag > 0.0) ? GCS_e_hPmag : LAS_e_hPmag )";delcols "GCS_hPmag GCS_e_hPmag LAS_hPmag LAS_e_hPmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize Petrosian J magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits cmd='addcol -units mag UK_jPmag "( (DXS_jPmag > 0.0) ? DXS_jPmag : LAS_jPmag )";addcol -units mag UK_e_jPmag "( (DXS_jPmag > 0.0) ? DXS_e_jPmag : LAS_e_jPmag )";delcols "DXS_jPmag DXS_e_jPmag LAS_jPmag LAS_e_jPmag";colmeta -name DXS_jPmag "UK_jPmag";colmeta -name DXS_e_jPmag "UK_e_jPmag";addcol -units mag UK_jPmag "( (GCS_jPmag > 0.0) ? GCS_jPmag : DXS_jPmag )";addcol -units mag UK_e_jPmag "( (GCS_jPmag > 0.0) ? GCS_e_jPmag : DXS_e_jPmag )";delcols "GCS_jPmag GCS_e_jPmag DXS_jPmag DXS_e_jPmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize Petrosian Y magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits cmd='addcol -units mag UK_yPmag "( (GCS_yPmag > 0.0) ? GCS_yPmag : LAS_yPmag )";addcol -units mag UK_e_yPmag "( (GCS_yPmag > 0.0) ? GCS_e_yPmag : LAS_e_yPmag )";delcols "GCS_yPmag GCS_e_yPmag LAS_yPmag LAS_e_yPmag"'
cp temp.fits temp_in.fits
echo 
echo Rename Petrosian Z magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=SUSS_UK_pet.fits cmd='colmeta -name UK_zPmag "GCS_zPmag";colmeta -name UK_e_zPmag "GCS_e_zPmag"'
