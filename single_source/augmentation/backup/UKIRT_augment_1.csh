echo 
echo Add petrosian photometry to CDS-upload-xmatch UKIDSS datasets  
date
echo
echo now DXS survey
topcat -stilts coneskymatch in=SUSS_DXS.fits ifmt=fits out=out.fits icmd='select (ra_x>0.0)' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/dxs9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pJmag e_pJmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
topcat -stilts tskymatch2 in1=SUSS_DXS.fits ifmt1=fits in2=out.fits out=out2.fits ra1='ra_x' dec1='dec_x' ra2='ra_x_2' dec2='dec_x_2' error=0.5 find=best1
topcat -stilts tpipe in=out2.fits out=SUSS_DXS_pet.fits cmd='delcols "ra_x_2 dec_x_2 GroupID GroupSize Separation"'
date
echo
echo now GCS survey
topcat -stilts coneskymatch in=SUSS_GCS.fits ifmt=fits out=out.fits icmd='select (ra_x>0.0)' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/gcs9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pZmag e_pZmag pYmag e_pYmag pJmag e_pJmag pHmag e_pHmag pKmag1 e_pKmag1";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
topcat -stilts tskymatch2 in1=SUSS_GCS.fits ifmt1=fits in2=out.fits out=out2.fits ra1='ra_x' dec1='dec_x' ra2='ra_x_2' dec2='dec_x_2' error=0.5 find=best1
topcat -stilts tpipe in=out2.fits out=SUSS_GCS_pet.fits cmd='delcols "ra_x_2 dec_x_2 GroupID GroupSize Separation"'
date
echo
echo now LAS survey
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out.fits icmd='select (ra_x>0.0)' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
topcat -stilts tskymatch2 in1=SUSS_LAS.fits ifmt1=fits in2=out.fits out=out2.fits ra1='ra_x' dec1='dec_x' ra2='ra_x_2' dec2='dec_x_2' error=0.5 find=best1
topcat -stilts tpipe in=out2.fits out=SUSS_LAS_pet.fits cmd='delcols "ra_x_2 dec_x_2 GroupID GroupSize Separation"'
date



echo
echo Split the LAS into chunks as it takes too long in one go
echo
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out1.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>0)&&(N_SUMMARY<=1000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out2.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>1000)&&(N_SUMMARY<=2000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out3.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>2000)&&(N_SUMMARY<=3000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out4.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>3000)&&(N_SUMMARY<=4000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out5.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>4000)&&(N_SUMMARY<=5000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out6.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>5000)&&(N_SUMMARY<=6000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out7.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>6000)&&(N_SUMMARY<=7000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out8.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>7000)&&(N_SUMMARY<=8000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out9.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>8000)&&(N_SUMMARY<=9000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out10.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>9000)&&(N_SUMMARY<=10000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
date
topcat -stilts coneskymatch in=SUSS_LAS.fits ifmt=fits out=out11.fits icmd='select ((ra_x>0.0)&&(N_SUMMARY>10000)&&(N_SUMMARY<=11000))' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/319/las9?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
topcat -stilts tcat in="out1.fits out2.fits out3.fits out4.fits out5.fits out6.fits out7.fits out8.fits out9.fits out10.fits out11.fits" out=out.fits 








date
echo
echo now GPS survey
topcat -stilts coneskymatch in=SUSS_GPS.fits ifmt=fits out=out.fits icmd='select (ra_x>340.0)' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://vizier.cds.unistra.fr/viz-bin/conesearch/II/316/gps6?' sr=0.0001388888 verb=3 parallel=5 erract=retry4 ocmd='\
keepcols "ra_x dec_x pYmag e_pYmag pJmag1 e_pJmag1 pHmag e_pHmag pKmag e_pKmag";\
colmeta -name ra_x_2 "ra_x";\
colmeta -name dec_x_2 "dec_x";\
'
topcat -stilts tskymatch2 in1=SUSS_GPS.fits ifmt1=fits in2=out.fits out=out2.fits ra1='ra_x' dec1='dec_x' ra2='ra_x_2' dec2='dec_x_2' error=0.5 find=best1
topcat -stilts tpipe in=out2.fits out=SUSS_GPS_pet.fits cmd='delcols "ra_x_2 dec_x_2 GroupID GroupSize Separation"'
date



















topcat -stilts tapskymatch in=SUSS_GCS.fits ifmt=fits out=out.fits icmd='select (ra_x>0.0)' find=each ra=ra_x dec=dec_x servicetype=cone serviceurl='http://wfaudata.roe.ac.uk/ukidssDR9-dsa/DirectCone?DSACAT=UKIDSS_DR9&DSATAB=gcsDetection&' sr=0.0001388888 verb=3 parallel=5


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


