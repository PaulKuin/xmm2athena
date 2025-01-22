#!/bin/csh
#  {obs_epoch} {gtemp} 
if ($2 == "") then
   echo This script from Paul is to add epch positions to subset high PM gaia dr3
   echo parameters:  obsepoch gaiatempfile gaiaoutputfile
endif
#
set OBSEPOCH=$1
set IN=$2
set OUT="ogtemp.tmp"
echo observed epoch $OBSEPOCH
echo input file $IN
# add epoch 2000 positions
java -jar /Users/kuin/bin/topcat-full.jar -stilts tpipe \
in=$IN ifmt=fits out=$OUT ofmt=fits \
cmd="addcol ep2000 epochProp(2000-2016.0,ra,dec,parallax,pmra,pmdec,radial_velocity);\
addcol obsepoch1 epochProp($OBSEPOCH-2016.0,ra,dec,parallax,pmra,pmdec,radial_velocity);\
addcol ra2000ep pick(ep2000,0)[0];\
addcol dec2000ep pick(ep2000,1)[0];\
addcol raobsep pick(obsepoch1,0)[0];\
addcol decobsep pick(obsepoch1,1)[0];\
delcols 'ra dec';\
addcol -before 1 ra raobsep;\
addcol -before 2 dec decobsep;\
delcols 'ep2000 obsepoch1 *_gspphot has* *rv* phot* random* astrometr* pseudo* ipd* ruw* l b ';\
addcol obsEpoch $OBSEPOCH "
#
mv $OUT $IN
#
#
# changes Aug 10, 2022: delcols ra,dec, and then copy *obsep into them; add obsEpoch  