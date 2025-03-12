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
set temp1=$WORK/temp1.fits
set temp2=$WORK/temp2.fits
set sdss=$WORK/sdss.fits
echo
# prep: match the singlesource2 to the SUSS_aug files (needed if the aux match used just 
# the columns IAUNAME, OBSID, SRCNUM, obs_epoch, ra2000Ep, dec2000Ep 
echo
set INprep = $WORK/SUSS_Gaia_aug_V5AB.fits
set OUTprep = $WORK/SUSS_Gaia_aug_V5ABc.fits
echo
echo     make Stars training set using $INprep
echo
topcat -stilts tpipe in=$INprep ifmt=fits out=$temp1 ofmt=fits cmd="addcol pm_error sqrt(e_pmRA*e_pmRA+e_pmDE*e_pmDE)" 
topcat -stilts tpipe in=$temp1 ifmt=fits out=stars_training_setIDs.fits  cmd='select (PM/pm_error>30)&&((UVW2_QUALITY_FLAG==0)|(UVM2_QUALITY_FLAG==0)|(UVW1_QUALITY_FLAG==0)|(U_QUALITY_FLAG==0)|(B_QUALITY_FLAG==0)|(V_QUALITY_FLAG==0)&&(RPlx>10))' cmd='keepcols IAUNAME'
echo
cp $INprep $temp1
echo add col class0
echo
topcat -stilts tmatch2 in1=$temp1 in2=stars_training_setIDs.fits ifmt1=fits ifmt2=fits matcher=exact values1="IAUNAME" values2="IAUNAME2" join=all1 find=best1 icmd2='colmeta -name IAUNAME2 "IAUNAME"' ocmd='addcol class0 (equals(IAUNAME,IAUNAME2)?0:NULL)' out=$OUTprep ofmt=fits
echo
echo find match SUSS to V/154/sdss16; retain spCL and Q parameters
echo
topcat -stilts cdsskymatch cdstable=V/154/sdss16 in=$INprep out=$sdss ofmt=fits find=best ra=ra2000Ep dec=dec2000Ep radius=3 ocmd='keepcols "IAUNAME spCl Q"'
echo
echo     make QSO training set 
echo      
topcat -stilts tpipe in=$sdss ifmt=fits out=qso_training_setIDs.fits ofmt=fits cmd='select (length(spCl)==3)&&(Q==3)' 
echo
echo  add col class1
echo
topcat -stilts tmatch2 in1=$OUTprep in2=qso_training_setIDs.fits ifmt1=fits ifmt2=fits matcher=exact values1="IAUNAME" values2="IAUNAME3" join=all1 find=best1 icmd1="delcols IAUNAME2" icmd2='colmeta -name IAUNAME3 "IAUNAME"' ocmd='addcol class1 (equals(IAUNAME,IAUNAME3)?1:NULL)' out=$temp2 ofmt=fits
cp $temp2 $OUTprep 
echo
echo     make galaxy training set
echo
topcat -stilts tpipe in=$sdss ifmt=fits out=galaxy_training_setIDs.fits cmd='select (length(spCL)==6)&&(Q==3)' 
echo
echo add col class2
echo
topcat -stilts tmatch2 in1=$OUTprep in2=galaxy_training_setIDs.fits ifmt1=fits ifmt2=fits matcher=exact values1="IAUNAME" values2="IAUNAME2" join=all1 find=best1 icmd1='delcols IAUNAME3' icmd2='colmeta -name IAUNAME2 "IAUNAME"' ocmd='addcol class2 (equals(IAUNAME,IAUNAME2)?2:NULL)' out=$temp2 ofmt=fits
cp $temp2 $OUTprep
echo
echo merge to make the column class and clean up
echo
topcat -stilts tpipe in=$temp2 out=$OUTprep ifmt=fits ofmt=fits find=each cmd='addcol class class0+class1+class2 ' cmd='delcols "spCl_* Q_* IAUNAME3 class0 class1 class2"'
echo
