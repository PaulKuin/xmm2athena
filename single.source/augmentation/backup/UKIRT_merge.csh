echo 
echo Combine UKIRT DXS and GCS data.
date
echo 
topcat -stilts tmatch2 in1=SUSS_GCS.fits ifmt1=fits in2=SUSS_DXS.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd1=' \
colmeta -name GCS_zmag "zAperMag3"; \
colmeta -name GCS_e_zmag "zAperMag3Err"; \
colmeta -name GCS_ymag "yAperMag3"; \
colmeta -name GCS_e_ymag "yAperMag3Err"; \
colmeta -name GCS_jmag "jAperMag3"; \
colmeta -name GCS_e_jmag "jAperMag3Err"; \
colmeta -name GCS_hmag "hAperMag3"; \
colmeta -name GCS_e_hmag "hAperMag3Err"; \
colmeta -name GCS_kmag "k_1AperMag3"; \
colmeta -name GCS_e_kmag "k_1AperMag3Err"; \
delcols "JName ra_x dec_x sourceID mode epoch mergedClass angDist" \
' icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; \
colmeta -name DXS_jmag "jAperMag3"; \
colmeta -name DXS_e_jmag "jAperMag3Err"; \
colmeta -name DXS_kmag "kAperMag3"; \
colmeta -name DXS_e_kmag "kAperMag3Err"; \
keepcols "SRCNUM2 OBSID2 DXS_jmag DXS_e_jmag DXS_kmag DXS_e_kmag"' \
ocmd='delcols "SRCNUM2 OBSID2 GroupID GroupSize Separation"'
cp temp.fits temp_in.fits
echo 
echo add GPS survey
date
echo 
topcat -stilts tmatch2 in1=temp_in.fits ifmt1=fits in2=SUSS_GPS.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; \
colmeta -name GPS_jmag "jAperMag3"; \
colmeta -name GPS_e_jmag "jAperMag3Err"; \
colmeta -name GPS_hmag "hAperMag3"; \
colmeta -name GPS_e_hmag "hAperMag3Err"; \
colmeta -name GPS_kmag "k_1AperMag3"; \
colmeta -name GPS_e_kmag "k_1AperMag3Err"; \
colmeta -name GP2_kmag "k_2AperMag3"; \
colmeta -name GP2_e_kmag "k_2AperMag3Err"; \
keepcols "SRCNUM2 OBSID2 GPS_jmag GPS_e_jmag GPS_hmag GPS_e_hmag GPS_kmag GPS_e_kmag GP2_kmag GP2_e_kmag"' \
ocmd='delcols "SRCNUM2 OBSID2 GroupID GroupSize Separation"'
cp temp.fits temp_in.fits
echo 
echo add LAS survey
date
echo 
topcat -stilts tmatch2 in1=temp_in.fits ifmt1=fits in2=SUSS_LAS.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; \
colmeta -name LAS_ymag "yAperMag3"; \
colmeta -name LAS_e_ymag "yAperMag3Err"; \
colmeta -name LAS_jmag "j_1AperMag3"; \
colmeta -name LAS_e_jmag "j_1AperMag3Err"; \
colmeta -name LAS_hmag "hAperMag3"; \
colmeta -name LAS_e_hmag "hAperMag3Err"; \
colmeta -name LAS_kmag "kAperMag3"; \
colmeta -name LAS_e_kmag "kAperMag3Err"; \
keepcols "SRCNUM2 OBSID2 LAS_ymag LAS_e_ymag LAS_jmag LAS_e_jmag LAS_hmag LAS_e_hmag LAS_kmag LAS_e_kmag"' \
ocmd='delcols "SRCNUM2 OBSID2 GroupID GroupSize Separation"'
cp temp.fits temp_in.fits
echo 
echo Prioritize K magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits \
cmd='addcol -units mag UK_kmag "( (GPS_kmag > 0.0) ? GPS_kmag : GP2_kmag )";\
addcol -units mag UK_e_kmag "( (GPS_kmag > 0.0) ? GPS_e_kmag : GP2_e_kmag )";\
delcols "GPS_kmag GPS_e_kmag GP2_kmag GP2_e_kmag";\
colmeta -name GPS_kmag "UK_kmag";\
colmeta -name GPS_e_kmag "UK_e_kmag";\
addcol -units mag UK_kmag "( (GPS_kmag > 0.0) ? GPS_kmag : LAS_kmag )";\
addcol -units mag UK_e_kmag "( (GPS_kmag > 0.0) ? GPS_e_kmag : LAS_e_kmag )";\
delcols "GPS_kmag GPS_e_kmag LAS_kmag LAS_e_kmag";\
colmeta -name GPS_kmag "UK_kmag";\
colmeta -name GPS_e_kmag "UK_e_kmag";\
addcol -units mag UK_kmag "( (DXS_kmag > 0.0) ? DXS_kmag : GPS_kmag )";\
addcol -units mag UK_e_kmag "( (DXS_kmag > 0.0) ? DXS_e_kmag : GPS_e_kmag )";\
delcols "DXS_kmag DXS_e_kmag GPS_kmag GPS_e_kmag";\
colmeta -name DXS_kmag "UK_kmag";\
colmeta -name DXS_e_kmag "UK_e_kmag";\
addcol -units mag UK_kmag "( (GCS_kmag > 0.0) ? GCS_kmag : DXS_kmag )";\
addcol -units mag UK_e_kmag "( (GCS_kmag > 0.0) ? GCS_e_kmag : DXS_e_kmag )";\
delcols "GCS_kmag GCS_e_kmag DXS_kmag DXS_e_kmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize H magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits \
cmd='\
addcol -units mag UK_hmag "( (GPS_hmag > 0.0) ? GPS_hmag : LAS_hmag )";\
addcol -units mag UK_e_hmag "( (GPS_hmag > 0.0) ? GPS_e_hmag : LAS_e_hmag )";\
delcols "GPS_hmag GPS_e_hmag LAS_hmag LAS_e_hmag";\
colmeta -name GPS_hmag "UK_hmag";\
colmeta -name GPS_e_hmag "UK_e_hmag";\
addcol -units mag UK_hmag "( (GCS_hmag > 0.0) ? GCS_hmag : GPS_hmag )";\
addcol -units mag UK_e_hmag "( (GCS_hmag > 0.0) ? GCS_e_hmag : GPS_e_hmag )";\
delcols "GCS_hmag GCS_e_hmag GPS_hmag GPS_e_hmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize J magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits \
cmd='\
addcol -units mag UK_jmag "( (GPS_jmag > 0.0) ? GPS_jmag : LAS_jmag )";\
addcol -units mag UK_e_jmag "( (GPS_jmag > 0.0) ? GPS_e_jmag : LAS_e_jmag )";\
delcols "GPS_jmag GPS_e_jmag LAS_jmag LAS_e_jmag";\
colmeta -name GPS_jmag "UK_jmag";\
colmeta -name GPS_e_jmag "UK_e_jmag";\
addcol -units mag UK_jmag "( (DXS_jmag > 0.0) ? DXS_jmag : GPS_jmag )";\
addcol -units mag UK_e_jmag "( (DXS_jmag > 0.0) ? DXS_e_jmag : GPS_e_jmag )";\
delcols "DXS_jmag DXS_e_jmag GPS_jmag GPS_e_jmag";\
colmeta -name DXS_jmag "UK_jmag";\
colmeta -name DXS_e_jmag "UK_e_jmag";\
addcol -units mag UK_jmag "( (GCS_jmag > 0.0) ? GCS_jmag : DXS_jmag )";\
addcol -units mag UK_e_jmag "( (GCS_jmag > 0.0) ? GCS_e_jmag : DXS_e_jmag )";\
delcols "GCS_jmag GCS_e_jmag DXS_jmag DXS_e_jmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize Y magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits \
cmd='\
addcol -units mag UK_ymag "( (GCS_ymag > 0.0) ? GCS_ymag : LAS_ymag )";\
addcol -units mag UK_e_ymag "( (GCS_ymag > 0.0) ? GCS_e_ymag : LAS_e_ymag )";\
delcols "GCS_ymag GCS_e_ymag LAS_ymag LAS_e_ymag"'
cp temp.fits temp_in.fits
echo 
echo Rename Z magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=SUSS_UK.fits \
cmd='\
colmeta -name UK_zmag "GCS_zmag";\
colmeta -name UK_e_zmag "GCS_e_zmag"'


