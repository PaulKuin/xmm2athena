#!/bin/csh 
#
# 
#  stilts_joined_mags 
#   - remove columns not needed for the classification in the input catalogue
#   - add mag columns from version input catalogue which was matched to SDSS13 
#   - add columns with weighted value the magnitude
#   - remove the instrumental magnitudes once used to compute the mean
#  scripts in ~/bin
#  data files in /Users/data/catalogs subdirectories: see below 
#  output data cover the SED of the UV-optical-IR for SUSS sources, with added Gaia PM and parallax
#
# NPMK 2023/ XMM2Athena WG7
#$HOME=/Users/kuin
set BIN=$HOME/bin
set WORK=`pwd`
# from Atrium SUSS single source catalogue augmented with PanStars,SkyMapper, UKIDS, VISTA and ALLWISE
set IN = /Users/data/catalogs/atrium/XMMOM_SUSS5.0_Sources_v0_1_aug.fits
# IN2 is an intermediate product by merging SUSS and SDSS16 in TopCat using CDS crossmatch
set IN2 = /Users/data/catalogs/suss_gaia_epic/XMMOM_SUSS5.0_Sources_xSDSSxGaia_v0.1.csv
# output file (both FITS and CSV)
set OUT1 = /Users/data/catalogs/suss_gaia_epic/XMMOM_SUSS5.0_Sources_v0_1_traininginput.fits
set OUT = /Users/data/catalogs/suss_gaia_epic/XMMOM_SUSS5.0_Sources_v0_1_traininginput.csv
set temp1=$WORK/temp1.fit
set temp2=$WORK/temp2.fit
#
echo input augmented SUSS v0.1 Jan 2023: $IN
echo input SUSS with Gaia DR3 plx and PM, merged with SDSS16 intermediate product: $IN2
echo output file for classification : $OUT
echo " "
echo "merging the two inputs; removing those columns not needed for output creation"
echo "i.e., removing from consideration SUSS+ RA_??? ?II *_FLUX* *_ST *pmra* *pmdec* *radial*"
echo "and RA_HMS DEC_DMS lII bII *_FLUX* *_ST *pmra* *pmdec* *radial* from SUSSxSDSS, then tmatch2 on IAUNAME"
java -jar $BIN/topcat-full.jar -stilts tmatch2 in1=$IN in2=$IN2 ifmt1=fits ifmt2=csv \
   icmd1="delcols  'RA_hms dec_dms OBSID SRCNUM *SRCDIST *SIGNIF *II *_FLUX* *_ST *pmra* *pmdec* *radial*' " \
   icmd2="delcols  'OBSID SRCNUM *_SRCDIST RA* DEC* POSERR lII bII *SIGNIF* *_FLUX* *_AB_MAG* *_FLAG* gaiadr3* *_VAR objID *sp* *radial* SDSS16'" \
   matcher=exact  values1=IAUNAME   values2=IAUNAME  \
   join=all1 find=best1 out=$temp2 ofmt=fits  ocmd='delcols IAUNAME_2'\
   fixcols=dups suffix1=  suffix2=_2 
 

echo "run merging weighted data from csh-scripts/stilts-joined* to merge bands"
#
# gmag (add PS, SM, SDSS)
$BIN/stilts_joined_mags3 $temp2 $temp1 gmag PS_gmag PS_e_gmag SM_gmag SM_e_gmag gmag e_gmag
echo gmag \(add PS, SM, SDSS\) done \n
exit
#
# gPKmag Kron and Petrosian (PS, SM)
$BIN/stilts_joined_mags2 $temp1 $temp2 gPKmag PS_gKmag PS_e_gKmag SM_gPmag SM_e_gPmag 
echo gPKmag Kron and Petrosian \(PS, SM\) done
echo " "
#
# umag (add OM, SM, SDSS) even though they have different band centres ...
$BIN/stilts_joined_mags3 $temp2 $temp1 umag SM_umag SM_e_gmag U_AB_MAG U_AB_MAG_ERR umag e_umag
echo umag \(add OM+SM\) done
echo " "
#
# uPmag = SM_uPmag (i.e., not renamed)
#
# vmag (add OM and SM)
$BIN/stilts_joined_mags2 $temp1 $temp2 vmag V_AB_MAG V_AB_MAG_ERR  SM_vmag SM_e_vmag 
echo vmag \(add OM and SM\) done
echo " "
#
# vPmag = SM_vPamg (i.e., not renamed)
#
# rmag (add PS, SM, SDSS)
$BIN/stilts_joined_mags3 $temp2 $temp1 rmag PS_rmag PS_e_rmag SM_rmag SM_e_rmag rmag e_rmag
echo rmag \(add PS, SM, SDSS\) done
#
# rPKmag (add PS Kron and SM Petrosian)
$BIN/stilts_joined_mags2 $temp1 $temp2 rPKmag PS_rKmag PS_e_rKmag SM_rPmag SM_e_rPmag 
echo rPKmag \(add PS Kron and SM Petrosian\) done
#
# imag (add PS, SM, SDSS)
$BIN/stilts_joined_mags3 $temp2 $temp1 imag PS_imag PS_e_imag SM_imag SM_e_imag imag e_imag
 imag \(add PS, SM, SDSS\) done
# 
# iPKmag (add PS Kron and SM Petrosian)
echo$BIN/stilts_joined_mags2 $temp1 $temp2 iPKmag PS_iKmag PS_e_iKmag SM_iPmag SM_e_iPmag
echo iPKmag \(add PS Kron and SM Petrosian\) done
#
# zmag (add PS, SM, UK, VI, SDSS)
$BIN/stilts_joined_mags5 $temp2 $temp1 zmag PS_zmag PS_e_zmag SM_zmag SM_e_zmag UK_zmag UK_e_zmag VI_zmag VI_e_zmag zmag e_zmag
echo zmag \(add PS, SM, UK, VI, SDSS\) done
echo " "
#
# zPKmag (add PS Kron and SM, VI Petrosian)
$BIN/stilts_joined_mags2 $temp1 $temp2 zPKmag PS_zKmag PS_e_zKmag SM_zPmag SM_e_zPmag VI_zPmag VI_e_zPmag 
echo zPKmag \(add PS Kron and SM, VI Petrosian\) done
#
# ymag (add PS, UK, VI)
$BIN/stilts_joined_mags3 $temp2 $temp1 ymag PS_ymag PS_e_ymag UK_ymag UK_e_ymag VI_ymag VI_e_ymag
echo ymag \(add PS, UK, VI\) done
#
# yPKmag (add PS Kron and SM, VI Petrosian)
$BIN/stilts_joined_mags3 $temp1 $temp2 yPKmag PS_yKmag PS_e_yKmag SM_yPmag SM_e_yPmag VI_yPmag VI_e_yPmag
echo yPKmag \(add PS Kron and SM, VI Petrosian\) done
#
# jmag (add UK, VI)
$BIN/stilts_joined_mags2 $temp2 $temp1 jmag  UK_jmag UK_e_jmag VI_jmag VI_e_jmag
echo jmag \(add UK, VI\) done
echo " "
#
# jPmag = UK_jPmag (i.e., not renamed)
#
# hmag (add UK, VI)
$BIN/stilts_joined_mags2 $temp1 $temp2 hmag  UK_hmag UK_e_hmag VI_hmag VI_e_hmag
echo hmag \(add UK, VI\) done
#
# hPmag = UK_hPmag (i.e., not renamed)
#
# kmag (add UK, VI)
$BIN/stilts_joined_mags2 $temp2 $temp1 kmag  UK_kmag UK_e_kmag VI_kmag VI_e_kmag
echo kmag \(add UK, VI\) done
echo "\n-------------\n "
#
# kPmag = UK_kPmag (i.e., not renamed)
#
# no joins for WISE mag's WI_W1mag,WI_W2mag,WI_W3mag,WI_W4mag
# no joins for UVOT uvw2, uvm2, uvw1, b
#
# Question: should we add old catalogue with B-magnitudes like USNO-B1 ? Probably not.
#
#   remove catalogue mag columns no longer needed (keep  uvw1, gmag, zmag, WI_w1mag)
#   define colour terms 
#   remove merged mag columns no longer needed for classification
#
java -jar $BIN/topcat-full.jar -stilts tpipe \
  in=$temp1 \
  out=$OUT1  \
  ifmt='fits' ofmt='fits' \ 
  cmd=" delcols 'PS_* SM_?mag SM_e_?mag UK_?mag UK_e_?mag VI_* '; \ 
  addcol uvw2_uvm2 UVW2_AB_MAG-UVM2_AB_MAG; \  
  addcol uvm2_uvw1 UVM2_AB_MAG-UVW1_AB_MAG; \  
  addcol uvw1_u    UVW1_AB_MAG-umag;\
  addcol b_v       B_AB_MAG-vmag; \
  addcol g_r       gmag-rmag; \
  addcol r_i       rmag-imag; \
  addcol i_z       imag-zmag; \
  addcol z-y       zmag-ymag; \
  addcol j_z       jmag-zmag; \
  addcol h_j       hmag-jmag; \
  addcol k_h       kmag-hmag; \
  addcol W1_z      WI_W1mag-zmag; \
  addcol W2_W1     WI_W2mag-WI_W1mag; \
  addcol W3_w2     WI_W3mag-WI_W2mag; \
  addcol W4_w3     WI_W4mag-WI_W3mag; \
  addcol gPK_uPmag gPKmag-SM_uPmag;  \
  addcol gPK_g     pPKmag-gmag;  \
  addcol rPK_gPmag rPKmag-gPKmag;  \
  addcol rPK_r     pPKmag-gmag;  \
  addcol iPK_rPmag iPKmag-rPKmag;  \
  addcol iPK_i     iPKmag-imag;  \
  addcol zPK_iPmag zPKmag-iPKmag;  \
  addcol zPK_z     zPKmag-zmag;  \
  delcols  imag jmag hmag kmag WI_e_* WI_W4mag WI_W3mag WI_W2mag "

echo FITS output done; CSV next  
java -jar $BIN/topcat-full.jar -stilts tpipe \
  in=$OUT1 \
  out=$OUT  \
  ifmt='fits' ofmt='csv' 
    
echo colours and normal-Petrosion/Kron done and useless magnitudes removed
echo colours are defined using underscores to separate bands
echo imag,jmag,hmag,kmag,WI_e_*,WI_W4mag,WI_W3mag,WI_W2mag removed
echo -----
# rm temp1.fit temp2.fit
echo Completed merging/cleaning/colours script