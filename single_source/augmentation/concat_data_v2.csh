#!/bin/csh
alias topcat 'java -jar ~/bin/topcat-full.jar '
echo 
echo Add PannSTARRS first.
date
echo 
topcat -stilts tmatch2 in1=singlesourcecat2.fits ifmt1=fits in2=SUSS_PS.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; \
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
date
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
cp temp.fits temp_in.fits
echo 
echo Add UKIDSS next.
date
echo 
topcat -stilts tmatch2 in1=temp_in.fits ifmt1=fits in2=SUSS_UK_pet.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd2='colmeta -name SRCNUM2 "SRCNUM"; \
 colmeta -name OBSID2 "OBSID"; keepcols "SRCNUM2 OBSID2 UK_zmag UK_e_zmag UK_ymag UK_e_ymag UK_jmag UK_e_jmag UK_hmag UK_e_hmag UK_kmag UK_e_kmag UK_zPmag UK_e_zPmag UK_yPmag UK_e_yPmag UK_jPmag UK_e_jPmag UK_hPmag UK_e_hPmag UK_kPmag UK_e_kPmag"' ocmd='delcols "SRCNUM2 OBSID2 Separation"'
cp temp.fits temp_in.fits
echo 
echo Add VISTA next.
date
echo 
topcat -stilts tmatch2 in1=temp_in.fits ifmt1=fits in2=SUSS_VI.fits out=temp.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; \
keepcols "SRCNUM2 OBSID2 VI_zmag VI_e_zmag VI_ymag VI_e_ymag VI_jmag VI_e_jmag VI_hmag VI_e_hmag VI_kmag VI_e_kmag"' ocmd='delcols "SRCNUM2 OBSID2 GroupID GroupSize Separation"'
cp temp.fits temp_in.fits
echo 
echo Finally, add WISE.
date
echo 
topcat -stilts tmatch2 in1=temp_in.fits ifmt1=fits in2=SUSS_WIS.fits out=SUSS_aug_V5.fits matcher=exact+exact values1='SRCNUM OBSID' values2='SRCNUM2 OBSID2' find=best1 join=all1 fixcols=none icmd2='colmeta -name SRCNUM2 "SRCNUM"; colmeta -name OBSID2 "OBSID"; \
colmeta -name WI_W1mag "W1mag"; \
colmeta -name WI_e_W1mag "e_W1mag"; \
colmeta -name WI_W2mag "W2mag"; \
colmeta -name WI_e_W2mag "e_W2mag"; \
colmeta -name WI_W3mag "W3mag"; \
colmeta -name WI_e_W3mag "e_W3mag"; \
colmeta -name WI_W4mag "W4mag"; \
colmeta -name WI_e_W4mag "e_W4mag"; \
keepcols "SRCNUM2 OBSID2 WI_W1mag WI_e_W1mag WI_W2mag WI_e_W2mag WI_W3mag WI_e_W3mag WI_W4mag WI_e_W4mag"' ocmd='delcols "SRCNUM2 OBSID2 GroupID GroupSize Separation"'


