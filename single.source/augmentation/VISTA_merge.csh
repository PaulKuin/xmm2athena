echo 
echo Combine VISTA Viking and VHS data.
date
echo 
topcat -stilts tmatch2 in1=SUSS_VIK.fits ifmt1=fits in2=SUSS_VHS.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd1=' \
colmeta -name VIK_zmag "Zpmag"; \
colmeta -name VIK_e_zmag "e_Zpmag"; \
colmeta -name VIK_ymag "Ypmag"; \
colmeta -name VIK_e_ymag "e_Ypmag"; \
colmeta -name VIK_jmag "Jpmag"; \
colmeta -name VIK_e_jmag "e_Jpmag"; \
colmeta -name VIK_hmag "Hpmag"; \
colmeta -name VIK_e_hmag "e_Hpmag"; \
colmeta -name VIK_kmag "Kspmag"; \
colmeta -name VIK_e_kmag "e_Kspmag"; \
colmeta -name VIK_class "Mclass"; \
delcols "SrcID RAdeg DEdeg name PriOrSec E(B-V) Zap3 e_Zap3 Zperrbits Yap3 e_Yap3 Yperrbits Jap3 e_Jap3 Jperrbits Hap3 e_Hap3 Hperrbits Ksap3 e_Ksap3 Ksperrbits angDist" \
' icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; \
colmeta -name VHS_ymag "Ypmag"; \
colmeta -name VHS_e_ymag "e_Ypmag"; \
colmeta -name VHS_jmag "Jpmag"; \
colmeta -name VHS_e_jmag "e_Jpmag"; \
colmeta -name VHS_hmag "Hpmag"; \
colmeta -name VHS_e_hmag "e_Hpmag"; \
colmeta -name VHS_kmag "Kspmag"; \
colmeta -name VHS_e_kmag "e_Kspmag"; \
colmeta -name VHS_class "Mclass"; \
keepcols "SRCNUM2 OBSID2 VHS_ymag VHS_e_ymag VHS_jmag VHS_e_jmag VHS_hmag VHS_e_hmag VHS_kmag VHS_e_kmag VHS_class"' \
ocmd='delcols "SRCNUM2 OBSID2 GroupID GroupSize Separation VIK_class VHS_class"'
cp temp.fits temp_in.fits
echo 
echo 
echo Rename Z magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits \
cmd='\
colmeta -name VI_zmag "VIK_zmag";\
colmeta -name VI_e_zmag "VIK_e_zmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize Y magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits \
cmd='\
addcol -units mag VI_ymag "( (VIK_ymag > 0.0) ? VIK_ymag : VHS_ymag )";\
addcol -units mag VI_e_ymag "( (VIK_ymag > 0.0) ? VIK_e_ymag : VHS_e_ymag )";\
delcols "VIK_ymag VIK_e_ymag VHS_ymag VHS_e_ymag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize J magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits \
cmd='\
addcol -units mag VI_jmag "( (VIK_jmag > 0.0) ? VIK_jmag : VHS_jmag )";\
addcol -units mag VI_e_jmag "( (VIK_jmag > 0.0) ? VIK_e_jmag : VHS_e_jmag )";\
delcols "VIK_jmag VIK_e_jmag VHS_jmag VHS_e_jmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize H magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits \
cmd='\
addcol -units mag VI_hmag "( (VIK_hmag > 0.0) ? VIK_hmag : VHS_hmag )";\
addcol -units mag VI_e_hmag "( (VIK_hmag > 0.0) ? VIK_e_hmag : VHS_e_hmag )";\
delcols "VIK_hmag VIK_e_hmag VHS_hmag VHS_e_hmag"'
cp temp.fits temp_in.fits
echo 
echo Prioritize K magnitudes.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=SUSS_VI.fits \
cmd='\
addcol -units mag VI_kmag "( (VIK_kmag > 0.0) ? VIK_kmag : VHS_kmag )";\
addcol -units mag VI_e_kmag "( (VIK_kmag > 0.0) ? VIK_e_kmag : VHS_e_kmag )";\
delcols "VIK_kmag VIK_e_kmag VHS_kmag VHS_e_kmag"'
