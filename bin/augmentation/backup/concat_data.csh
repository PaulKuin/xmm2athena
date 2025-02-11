echo 
echo Add PannSTARRS first.
echo 
topcat -stilts tmatch2 in1=XMM-OM-SUSS5.0.fits.gz ifmt1=fits in2=SUSS_PS.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; \
colmeta -name PS_gmag "gmag"; \
colmeta -name PS_e_gmag "e_gmag"; \
colmeta -name PS_gKmag "gKmag"; \
colmeta -name PS_e_gKmag "e_gKmag"; \
colmeta -name PS_rmag "rmag"; \
colmeta -name PS_e_rmag "e_rmag"; \
colmeta -name PS_rKmag "rKmag"; \
colmeta -name PS_e_rKmag "e_rKmag"; \
colmeta -name PS_imag "imag"; \
colmeta -name PS_e_imag "e_imag"; \
colmeta -name PS_iKmag "iKmag"; \
colmeta -name PS_e_iKmag "e_iKmag"; \
colmeta -name PS_zmag "zmag"; \
colmeta -name PS_e_zmag "e_zmag"; \
colmeta -name PS_zKmag "zKmag"; \
colmeta -name PS_e_zKmag "e_zKmag"; \
colmeta -name PS_ymag "ymag"; \
colmeta -name PS_e_ymag "e_ymag"; \
colmeta -name PS_yKmag "yKmag"; \
colmeta -name PS_e_yKmag "e_yKmag"; \
keepcols "SRCNUM2 OBSID2 PS_gmag PS_e_gmag PS_gKmag PS_e_gKmag PS_rmag PS_e_rmag PS_rKmag PS_e_rKmag PS_imag PS_e_imag PS_iKmag  PS_e_iKmag PS_zmag PS_e_zmag PS_zKmag PS_e_zKmag PS_ymag PS_e_ymag PS_yKmag PS_e_yKmag"' ocmd='delcols "SRCNUM2 OBSID2 GroupID GroupSize Separation"'
cp temp.fits temp_in.fits
echo 
echo Add SkyMapper next.
echo 
topcat -stilts tmatch2 in1=temp_in.fits ifmt1=fits in2=SUSS_SM.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; \
colmeta -name SM_stellarity "ClassStar"; \
colmeta -name SM_umag "uPSF"; \
colmeta -name SM_e_umag "e_uPSF"; \
colmeta -name SM_uPmag "uPetro"; \
colmeta -name SM_e_uPmag "e_uPetro"; \
colmeta -name SM_vmag "vPSF"; \
colmeta -name SM_e_vmag "e_vPSF"; \
colmeta -name SM_vPmag "vPetro"; \
colmeta -name SM_e_vPmag "e_vPetro"; \
colmeta -name SM_gmag "gPSF"; \
colmeta -name SM_e_gmag "e_gPSF"; \
colmeta -name SM_gPmag "gPetro"; \
colmeta -name SM_e_gPmag "e_gPetro"; \
colmeta -name SM_rmag "rPSF"; \
colmeta -name SM_e_rmag "e_rPSF"; \
colmeta -name SM_rPmag "rPetro"; \
colmeta -name SM_e_rPmag "e_rPetro"; \
colmeta -name SM_imag "iPSF"; \
colmeta -name SM_e_imag "e_iPSF"; \
colmeta -name SM_iPmag "iPetro"; \
colmeta -name SM_e_iPmag "e_iPetro"; \
colmeta -name SM_zmag "zPSF"; \
colmeta -name SM_e_zmag "e_zPSF"; \
colmeta -name SM_zPmag "zPetro"; \
colmeta -name SM_e_zPmag "e_zPetro"; \
keepcols "SRCNUM2 OBSID2 SM_umag SM_e_umag SM_uPmag SM_e_uPmag SM_vmag SM_e_vmag SM_vPmag SM_e_vPmag SM_gmag SM_e_gmag SM_gPmag SM_e_gPmag SM_rmag SM_e_rmag SM_rPmag SM_e_rPmag SM_imag SM_e_imag SM_iPmag  SM_e_iPmag SM_zmag SM_e_zmag SM_zPmag SM_e_zPmag"' ocmd='delcols "SRCNUM2 OBSID2 GroupID GroupSize Separation"'
