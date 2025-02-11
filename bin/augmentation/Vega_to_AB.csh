echo 
echo Convert UKIDSS magnitudes to AB.
date
echo 
topcat -stilts tpipe in=SUSS_aug_V3.fits.gz out=temp.fits \
cmd='\
addcol -units ABmag AB_UK_zmag "(UK_zmag + 0.528)";\
addcol -units ABmag AB_UK_e_zmag "UK_e_zmag";\
delcols "UK_zmag UK_e_zmag";\
colmeta -name UK_zmag "AB_UK_zmag";\
colmeta -name UK_e_zmag "AB_UK_e_zmag";\
addcol -units ABmag AB_UK_ymag "(UK_ymag + 0.634)";\
addcol -units ABmag AB_UK_e_ymag "UK_e_ymag";\
delcols "UK_ymag UK_e_ymag";\
colmeta -name UK_ymag "AB_UK_ymag";\
colmeta -name UK_e_ymag "AB_UK_e_ymag";\
addcol -units ABmag AB_UK_jmag "(UK_jmag + 0.938)";\
addcol -units ABmag AB_UK_e_jmag "UK_e_jmag";\
delcols "UK_jmag UK_e_jmag";\
colmeta -name UK_jmag "AB_UK_jmag";\
colmeta -name UK_e_jmag "AB_UK_e_jmag";\
addcol -units ABmag AB_UK_hmag "(UK_hmag + 1.379)";\
addcol -units ABmag AB_UK_e_hmag "UK_e_hmag";\
delcols "UK_hmag UK_e_hmag";\
colmeta -name UK_hmag "AB_UK_hmag";\
colmeta -name UK_e_hmag "AB_UK_e_hmag";\
addcol -units ABmag AB_UK_kmag "(UK_kmag + 1.900)";\
addcol -units ABmag AB_UK_e_kmag "UK_e_kmag";\
delcols "UK_kmag UK_e_kmag";\
colmeta -name UK_kmag "AB_UK_kmag";\
colmeta -name UK_e_kmag "AB_UK_e_kmag";\
addcol -units ABmag AB_UK_zPmag "(UK_zPmag + 0.528)";\
addcol -units ABmag AB_UK_e_zPmag "UK_e_zPmag";\
delcols "UK_zPmag UK_e_zPmag";\
colmeta -name UK_zPmag "AB_UK_zPmag";\
colmeta -name UK_e_zPmag "AB_UK_e_zPmag";\
addcol -units ABmag AB_UK_yPmag "(UK_yPmag + 0.634)";\
addcol -units ABmag AB_UK_e_yPmag "UK_e_yPmag";\
delcols "UK_yPmag UK_e_yPmag";\
colmeta -name UK_yPmag "AB_UK_yPmag";\
colmeta -name UK_e_yPmag "AB_UK_e_yPmag";\
addcol -units ABmag AB_UK_jPmag "(UK_jPmag + 0.938)";\
addcol -units ABmag AB_UK_e_jPmag "UK_e_jPmag";\
delcols "UK_jPmag UK_e_jPmag";\
colmeta -name UK_jPmag "AB_UK_jPmag";\
colmeta -name UK_e_jPmag "AB_UK_e_jPmag";\
addcol -units ABmag AB_UK_hPmag "(UK_hPmag + 1.379)";\
addcol -units ABmag AB_UK_e_hPmag "UK_e_hPmag";\
delcols "UK_hPmag UK_e_hPmag";\
colmeta -name UK_hPmag "AB_UK_hPmag";\
colmeta -name UK_e_hPmag "AB_UK_e_hPmag";\
addcol -units ABmag AB_UK_kPmag "(UK_kPmag + 1.900)";\
addcol -units ABmag AB_UK_e_kPmag "UK_e_kPmag";\
delcols "UK_kPmag UK_e_kPmag";\
colmeta -name UK_kPmag "AB_UK_kPmag";\
colmeta -name UK_e_kPmag "AB_UK_e_kPmag";\
'
cp temp.fits temp_in.fits
echo 
echo Convert VISTA mags to AB next.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits \
cmd='\
addcol -units ABmag AB_VI_zmag "(VI_zmag + 0.502)";\
addcol -units ABmag AB_VI_e_zmag "VI_e_zmag";\
delcols "VI_zmag VI_e_zmag";\
colmeta -name VI_zmag "AB_VI_zmag";\
colmeta -name VI_e_zmag "AB_VI_e_zmag";\
addcol -units ABmag AB_VI_ymag "(VI_ymag + 0.600)";\
addcol -units ABmag AB_VI_e_ymag "VI_e_ymag";\
delcols "VI_ymag VI_e_ymag";\
colmeta -name VI_ymag "AB_VI_ymag";\
colmeta -name VI_e_ymag "AB_VI_e_ymag";\
addcol -units ABmag AB_VI_jmag "(VI_jmag + 0.916)";\
addcol -units ABmag AB_VI_e_jmag "VI_e_jmag";\
delcols "VI_jmag VI_e_jmag";\
colmeta -name VI_jmag "AB_VI_jmag";\
colmeta -name VI_e_jmag "AB_VI_e_jmag";\
addcol -units ABmag AB_VI_hmag "(VI_hmag + 1.366)";\
addcol -units ABmag AB_VI_e_hmag "VI_e_hmag";\
delcols "VI_hmag VI_e_hmag";\
colmeta -name VI_hmag "AB_VI_hmag";\
colmeta -name VI_e_hmag "AB_VI_e_hmag";\
addcol -units ABmag AB_VI_kmag "(VI_kmag + 1.827)";\
addcol -units ABmag AB_VI_e_kmag "VI_e_kmag";\
delcols "VI_kmag VI_e_kmag";\
colmeta -name VI_kmag "AB_VI_kmag";\
colmeta -name VI_e_kmag "AB_VI_e_kmag";\
'
cp temp.fits temp_in.fits
echo 
echo Convert WISE mags to AB next.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=temp.fits \
cmd='\
addcol -units ABmag AB_WI_W1mag "(WI_W1mag + 2.699)";\
addcol -units ABmag AB_WI_e_W1mag "WI_e_W1mag";\
delcols "WI_W1mag WI_e_W1mag";\
colmeta -name WI_W1mag "AB_WI_W1mag";\
colmeta -name WI_e_W1mag "AB_WI_e_W1mag";\
addcol -units ABmag AB_WI_W2mag "(WI_W2mag + 3.339)";\
addcol -units ABmag AB_WI_e_W2mag "WI_e_W2mag";\
delcols "WI_W2mag WI_e_W2mag";\
colmeta -name WI_W2mag "AB_WI_W2mag";\
colmeta -name WI_e_W2mag "AB_WI_e_W2mag";\
addcol -units ABmag AB_WI_W3mag "(WI_W3mag + 5.174)";\
addcol -units ABmag AB_WI_e_W3mag "WI_e_W3mag";\
delcols "WI_W3mag WI_e_W3mag";\
colmeta -name WI_W3mag "AB_WI_W3mag";\
colmeta -name WI_e_W3mag "AB_WI_e_W3mag";\
addcol -units ABmag AB_WI_W4mag "(WI_W4mag + 6.620)";\
addcol -units ABmag AB_WI_e_W4mag "WI_e_W4mag";\
delcols "WI_W4mag WI_e_W4mag";\
colmeta -name WI_W4mag "AB_WI_W4mag";\
colmeta -name WI_e_W4mag "AB_WI_e_W4mag";\
'
cp temp.fits temp_in.fits
echo 
echo Now change the units of the SM and PS magnitudes to AB mag too.
date
echo 
topcat -stilts tpipe in=temp_in.fits out=SUSS_aug_V4.fits \
cmd='\
colmeta -units ABmag "PS_gmag";\
colmeta -units ABmag "PS_e_gmag";\
colmeta -units ABmag "PS_gKmag";\
colmeta -units ABmag "PS_e_gKmag";\
colmeta -units ABmag "PS_rmag";\
colmeta -units ABmag "PS_e_rmag";\
colmeta -units ABmag "PS_rKmag";\
colmeta -units ABmag "PS_e_rKmag";\
colmeta -units ABmag "PS_imag";\
colmeta -units ABmag "PS_e_imag";\
colmeta -units ABmag "PS_iKmag";\
colmeta -units ABmag "PS_e_iKmag";\
colmeta -units ABmag "PS_zmag";\
colmeta -units ABmag "PS_e_zmag";\
colmeta -units ABmag "PS_zKmag";\
colmeta -units ABmag "PS_e_zKmag";\
colmeta -units ABmag "PS_ymag";\
colmeta -units ABmag "PS_e_ymag";\
colmeta -units ABmag "PS_yKmag";\
colmeta -units ABmag "PS_e_yKmag";\
colmeta -units ABmag "SM_umag";\
colmeta -units ABmag "SM_e_umag";\
colmeta -units ABmag "SM_uPmag";\
colmeta -units ABmag "SM_e_uPmag";\
colmeta -units ABmag "SM_vmag";\
colmeta -units ABmag "SM_e_vmag";\
colmeta -units ABmag "SM_vPmag";\
colmeta -units ABmag "SM_e_vPmag";\
colmeta -units ABmag "SM_gmag";\
colmeta -units ABmag "SM_e_gmag";\
colmeta -units ABmag "SM_gPmag";\
colmeta -units ABmag "SM_e_gPmag";\
colmeta -units ABmag "SM_rmag";\
colmeta -units ABmag "SM_e_rmag";\
colmeta -units ABmag "SM_rPmag";\
colmeta -units ABmag "SM_e_rPmag";\
colmeta -units ABmag "SM_imag";\
colmeta -units ABmag "SM_e_imag";\
colmeta -units ABmag "SM_iPmag";\
colmeta -units ABmag "SM_e_iPmag";\
colmeta -units ABmag "SM_zmag";\
colmeta -units ABmag "SM_e_zmag";\
colmeta -units ABmag "SM_zPmag";\
colmeta -units ABmag "SM_e_zPmag";\
'
gzip SUSS_aug_V4.fits
