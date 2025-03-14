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
# NPMK 2023/ XMM2Athena WG7/
#  revised 2025-03 for use in XMM5 further SUSS processing (classifications) 
#
alias topcat 'java -jar ~/bin/topcat-full.jar'
#
#$HOME=/Users/kuin
set BIN=$HOME/bin
set WORK=`pwd`
#
set temp1=$WORK/temp1.fit
set temp2=$WORK/temp2.fit
echo
# prep: match the singlesource2 to the SUSS_aug files (needed if the aux match used just 
# the columns IAUNAME, OBSID, SRCNUM, obs_epoch, ra2000Ep, dec2000Ep 
echo 
#   fix the following after testing rest:
set INprep1 = $WORK/sussxgaiadr3_ep2000_singlerecs_stg2.fits
set INprep2 = $WORK/SUSS_aug_V5AB.fits
set OUTprep = $WORK/SUSS_Gaia_aug_V5AB.fits
java -jar $BIN/topcat-full.jar -stilts tmatch2 in1=$INprep1 in2=$IN2prep2 ifmt1=fits ifmt2=fits \
   icmd1="delcols  'RA_hms dec_dms *_FLUX*  ' " \
   icmd2="delcols  'OBSID SRCNUM obs_epoch ra2000Ep dec2000Ep'" \
   matcher=exact  values1="IAUNAME"   values2="IAUNAME"  \
   join=all1 find=best1 out=$temp1 ofmt=fits  ocmd='delcols "IAUNAME_2 GroupID GroupSize"'\
   colmeta -name IAUNAME "IAUNAME_1";\
   colmeta -name PM "PM_cds";\
   colmeta -name pmRA "pmRA_cds";\
   fixcols=dups suffix1=  suffix2=_2 
cp $temp1 $OUTprep   
echo
echo add class column for the training set ... add_class.csh calls add_class.py
add_class.csh
echo
#   end fix ---
#
# SUSS single source catalogue with Gaia DR3 match, augmented with PanStars,SkyMapper, UKIDS, VISTA and ALLWISE
echo
set IN=$WORK/SUSS_Gaia_aug_V5ABc.fits
echo
echo input SUSS + GaiaDR3 + augment March 2025: $IN
#
# output file (both FITS and CSV)
set OUT1=$HOME/SUSS_traininginput.fits
set OUT=$HOME/SUSS_traininginput.fits
#
echo output file for classification $OUT
echo 
echo simplifying by removing some columns for the classification file 
echo i.e., columns OBSID SRCNUM "SIGNIF*" "*FLUX*" "*_ST"
echo 
java -jar $BIN/topcat-full.jar -stilts tpipe in=$IN ifmt=fits out=$temp2 ofmt=fits cmd='delcols  "OBSID SRCNUM *SIGNIF *_ST";colmeta -name vPmag "SM_vPmag";colmeta -name jPmag "UK_jPmag";colmeta -name uPmag "SM_uPmag";colmeta -name hPmag "UK_hPmag";colmeta -name kPmag "UK_kPmag";'  
echo
echo run merging weighted data from csh-scripts/stilts-joined_mags.py to merge bands
echo 
#java -jar ~/bin/topcat-full.jar -stilts cdsskymatch cdstable=V/154/sdss16 find=each in=$temp2 ifmt=fits  ra=ra2000Ep dec=dec2000Ep radius=3 out=$IN /
# ocmd='delcols "objID *_ICRS errHalf* mode class_cds clean *zsp spCl subCl pmRA_cds pmDE_cds e_pmRA_cds e_pmDE_cds sig* zph e_zph <zph> Q SDSS16 Sp-ID MJD angDist" ' 

#
python stilts_joined_mags.py 
exit
#
#   remove catalogue mag columns no longer needed (keep  uvw1, gmag, zmag, WI_w1mag)
#   define colour terms for classification 
#   remove mag columns no longer needed for classification ? no
#
# see also email Mat Page Jan 31 2024 on gaia_extended
#    and April 19, 2024 on Gaia_G_Wise_W1 definitions
#
echo 
echo removing merged input columns
echo
java -jar $BIN/topcat-full.jar -stilts tpipe in=$temp1 out=$OUT1 ifmt='fits' ofmt='fits' cmd="\
delcols 'PS_* SM_?mag SM_e_?mag UK_?mag UK_e_?mag VI_*';\
addcol uvw2_uvw1 UVW2_AB_MAG-UVW1_AB_MAG;addcol uvm2_uvw1 UVM2_AB_MAG-UVW1_AB_MAG;\
addcol uvw1_u UVW1_AB_MAG-umag;addcol OMu-b U_AB_MAG-B_AB_MAG;addcol b_v bmag-vmag;\
addcol uvw1_Gmag UVW1_AB_MAG-Gmag;addcol u_r umag-rmag;addcol k_W1 kmag-WI_W1;\
addcol W2_W1 WI_W2mag-WI_W1mag;
addcol gaia_extended '((BII>10.0)||(BII<-10.0))&&(gmag_M<19.0)&&!(ra2000Ep<500.0)?(20.0-gmag_M):((gmag_M>0.0)?(Gmag-gmag_M):(Gmag-V_AB_MAG))';\
addcol Gaia_G_WISE_W1 '((BII>10.0)||(BII<-10.0))&&(WI_W1mag<16.0)&&!(ra2000Ep<500.0)?(20.0-WI_W1mag):(Gmag<15.6)&&!(WI_W1mag>0.0)?(Gmag-17.1):(Gmag-WI_W1mag)';\
colmeta -name BP_RP 'BP-RP';delcols 'imag jmag hmag kmag WI_e_* WI_W4mag WI_W3mag WI_W2mag';"
echo 
echo FITS output done; CSV next  
java -jar $BIN/topcat-full.jar -stilts tpipe \
  in=$OUT1 \
  out=$OUT  \
  ifmt='fits' ofmt='csv' 
    
echo colours and normal-Petrosion/Kron done and survey individual magnitudes removed
echo colours are defined using underscores to separate bands
echo imag,jmag,hmag,kmag,WI_e_*,WI_W4mag,WI_W3mag,WI_W2mag removed
echo -----
# rm temp1.fit temp2.fit
echo Completed merging/cleaning/colours script